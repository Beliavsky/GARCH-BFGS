! Fit joint overnight and intraday MCS-GARCH volatility models.
!
! The fitted likelihood includes close-to-open overnight returns and
! regular-session intraday close-to-close returns.  Overnight returns have an
! estimated variance weight, while intraday returns use a diurnal curve and
! MCS-GARCH intraday dynamics.

module fit_mcsgarch_on_intraday_mod
    use kind_mod, only: dp
    use math_const_mod, only: log_sqrt_2pi, pi
    use date_mod, only: date_label, yyyymmdd
    use market_data_mod, only: ohlcv_series_t, read_intraday_prices_file, filter_intraday_session, intraday_bin_ids
    use garch_mcsgarch_mod, only: mcsgarch_params_t, mcsgarch_fit_result_t, estimate_diurnal_variance, &
                                  mcsgarch_filter, mcsgarch_nagarch_filter, mcsgarch_gjr_filter
    use distributions_mod, only: pdf_fs_skewt
    use bfgs_mod, only: bfgs_minimize
    use stats_mod, only: mean, variance
    use input_files_mod, only: collect_input_filenames, MAX_PATH_LEN
    implicit none
    private

    character(len=*), parameter :: file_pattern = "c:\python\databento\data_1min\*.bin" ! "c:\python\intraday_prices\*.csv"
    real(dp), parameter :: min_var = 1.0e-12_dp
    real(dp), parameter :: min_pdf = 1.0e-300_dp
    real(dp), parameter :: persist_max = 0.999_dp
    logical, parameter :: smooth_diurnal_curve = .true.
    integer, parameter :: diurnal_smooth_half_width = 2
    integer, parameter :: model_sym = 1
    integer, parameter :: model_nagarch = 2
    integer, parameter :: model_gjr = 3
    integer, parameter :: dist_normal = 1
    integer, parameter :: dist_t = 2
    integer, parameter :: dist_fs_skewt = 3
    integer, parameter :: lambda_zero = 1
    integer, parameter :: lambda_est = 2
    integer, parameter :: lambda_one = 3
    character(len=12), parameter :: model_names(*) = [character(len=12) :: &
        "MCSGARCH", "MCSNAGARCH", "MCSGJRGARCH"]
    character(len=8), parameter :: dist_names(*) = [character(len=8) :: &
        "NORMAL", "T", "FS_SKEWT"]
    character(len=6), parameter :: lambda_names(*) = [character(len=6) :: &
        "ON0", "ONEST", "ON1"]

    type :: on_intraday_fit_t
        type(mcsgarch_fit_result_t) :: fit
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
    integer, save :: obj_model = model_sym
    integer, save :: obj_dist = dist_normal
    integer, save :: obj_lambda_mode = lambda_one

    public :: run_fit_mcsgarch_on_intraday

contains

    ! Fit all configured models to each input file obtained from command-line, glob, or default.
    subroutine run_fit_mcsgarch_on_intraday()
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
    end subroutine run_fit_mcsgarch_on_intraday

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
        call read_intraday_prices_file(trim(filename), bars)
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
        real(dp), allocatable :: p(:), grad(:), initial_intra_scale(:), fitted_intra_scale(:)
        real(dp) :: t0, t1
        integer :: model, dist, np, niter, nobs
        logical :: converged

        model = model_id(model_name)
        dist = dist_id(dist_name)
        np = model_param_count(model) + 1 + lambda_param_count(model, lambda_mode) + dist_param_count(dist)
        allocate(p(np), grad(np), initial_intra_scale(size(intra_ret)), fitted_intra_scale(size(intra_ret)))
        call intraday_scale_from_overnight(intra_base, intra_on_ret, model, initial_lambda_on(lambda_mode), &
                                           0.0_dp, initial_gamma_on(model, lambda_mode), initial_intra_scale)
        call estimate_diurnal_variance(intra_ret, initial_intra_scale, bin_id, diurnal_var, &
                                       smooth_diurnal_curve, diurnal_smooth_half_width)
        call set_objective_data(intra_ret, intra_base, intra_on_ret, diurnal_var, on_ret, on_scale, &
                                model, dist, lambda_mode)
        call pack_joint_params(initial_params(intra_ret, initial_intra_scale, diurnal_var), &
                               initial_overnight_weight(on_ret, on_scale), initial_lambda_on(lambda_mode), &
                               0.0_dp, initial_gamma_on(model, lambda_mode), model, dist, lambda_mode, p)

        call cpu_time(t0)
        call bfgs_minimize(joint_obj, p, np, max_iter, gtol, result%fit%nll, niter, converged)
        call cpu_time(t1)
        result%fit_sec = t1 - t0
        call unpack_joint_params(p, model, dist, result%fit%params, result%overnight_weight, &
                                 result%lambda_on, result%theta_on, result%gamma_on, result%fit%nu, result%fit%xi, &
                                 lambda_mode)
        call joint_obj(p, np, result%fit%nll, grad)
        call intraday_scale_from_overnight(intra_base, intra_on_ret, model, result%lambda_on, result%theta_on, &
                                           result%gamma_on, fitted_intra_scale)
        call filter_by_model(model, intra_ret, fitted_intra_scale, diurnal_var, result%fit%params, q)

        nobs = size(intra_ret) + size(on_ret)
        result%fit%loglik = -result%fit%nll * real(nobs, dp)
        result%fit%persist = persist(result%fit%params, model)
        result%fit%niter = niter
        result%fit%converged = converged
        result%aic = -2.0_dp*result%fit%loglik + 2.0_dp*real(model_param_count(model) + 1 + &
                     lambda_param_count(model, lambda_mode) + &
                     dist_param_count(dist) + count_populated_bins(bin_id), dp)
        result%bic = -2.0_dp*result%fit%loglik + log(real(nobs, dp)) * real(model_param_count(model) + 1 + &
                     lambda_param_count(model, lambda_mode) + &
                     dist_param_count(dist) + count_populated_bins(bin_id), dp)
        deallocate(p, grad, initial_intra_scale, fitted_intra_scale)
    end subroutine fit_one_joint_model

    ! Store data used by the BFGS objective.
    subroutine set_objective_data(intra_ret, intra_base, intra_on_ret, diurnal_var, on_ret, on_scale, &
                                  model, dist, lambda_mode)
        real(dp), intent(in) :: intra_ret(:), intra_base(:), intra_on_ret(:), diurnal_var(:), on_ret(:), on_scale(:)
        integer, intent(in) :: model, dist, lambda_mode

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
        obj_model = model
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
        allocate(pp(np), pm(np))
        do j = 1, np
            step = 1.0e-5_dp * max(1.0_dp, abs(p(j)))
            pp = p
            pm = p
            pp(j) = pp(j) + step
            pm(j) = pm(j) - step
            fp = joint_nll_from_p(pp)
            fm = joint_nll_from_p(pm)
            g(j) = (fp - fm) / (2.0_dp * step)
        end do
        deallocate(pp, pm)
    end subroutine joint_obj

    ! Return average negative log likelihood for packed parameters.
    real(dp) function joint_nll_from_p(p)
        real(dp), intent(in) :: p(:)
        type(mcsgarch_params_t) :: params
        real(dp), allocatable :: q(:), intra_scale(:)
        real(dp) :: on_weight, lambda_on, theta_on, gamma_on, nu, xi, loss, h
        integer :: i

        call unpack_joint_params(p, obj_model, obj_dist, params, on_weight, lambda_on, theta_on, gamma_on, nu, xi, &
                                 obj_lambda_mode)
        if (.not. params_valid(params, obj_model) .or. on_weight <= 0.0_dp .or. &
            lambda_on < 0.0_dp .or. gamma_on < 0.0_dp) then
            joint_nll_from_p = huge(1.0_dp) / 10.0_dp
            return
        end if
        allocate(q(size(obj_intra_ret)), intra_scale(size(obj_intra_ret)))
        call intraday_scale_from_overnight(obj_intra_base, obj_intra_on_ret, obj_model, lambda_on, theta_on, &
                                           gamma_on, intra_scale)
        call filter_by_model(obj_model, obj_intra_ret, intra_scale, obj_diurnal_var, params, q)
        loss = 0.0_dp
        do i = 1, size(obj_on_ret)
            h = max(obj_on_scale(i) * on_weight, min_var)
            loss = loss + innovation_nll(obj_on_ret(i), h, obj_dist, nu, xi)
        end do
        do i = 1, size(obj_intra_ret)
            h = max(intra_scale(i) * obj_diurnal_var(i) * q(i), min_var)
            loss = loss + innovation_nll(obj_intra_ret(i), h, obj_dist, nu, xi)
        end do
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

        print '(A)', "Joint overnight/intraday MCS-GARCH fits"
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
        print '(A12,1X,A8,1X,A6,1X,A12,11(1X,A10),1X,A8,1X,A5,1X,A13,1X,A9)', &
              "Model", "Dist", "Lambda", "omega", "alpha", "gamma", "beta", "theta", "on_wt", "lambda_on", "nu", "xi", &
              "theta_on", "gamma_on", "persist", "iter", "conv", "logL", "fit_sec"
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

        print '(A12,1X,A8,1X,A6,1X,ES12.4,11(1X,F10.4),1X,I8,5X,L1,1X,F13.3,1X,F9.3)', &
              model_name, dist_name, lambda_name, result%fit%params%omega, result%fit%params%alpha, &
              result%fit%params%gamma, result%fit%params%beta, result%fit%params%theta, &
              result%overnight_weight, result%lambda_on, result%fit%nu, result%fit%xi, result%theta_on, &
              result%gamma_on, result%fit%persist, &
              result%fit%niter, result%fit%converged, result%fit%loglik, result%fit_sec
    end subroutine print_fit_row

    ! Print the information criteria comparison table.
    subroutine print_ic_table(results)
        type(on_intraday_fit_t), intent(in) :: results(:)
        integer :: i

        print '(A)', ""
        print '(A)', "Joint overnight/intraday model comparison"
        print '(A)', "----------------------------------------------------------------------------------------------"
        print '(A18,1X,A8,1X,A6,1X,A10,1X,A10,1X,A10,1X,A14,1X,A14,1X,A14,1X,A8,1X,A8,1X,A8)', &
              "Model", "Dist", "Lambda", "lambda_on", "theta_on", "gamma_on", "logL", "AIC", "BIC", "conv", "iter", &
              "fit_sec"
        print '(A)', "----------------------------------------------------------------------------------------------"
        do i = 1, size(results)
            print '(A18,1X,A8,1X,A6,1X,F10.4,1X,F10.4,1X,F10.4,1X,F14.3,1X,F14.3,1X,F14.3,3X,L1,1X,I8,1X,F8.3)', &
                  model_table_label(trim(model_names(model_index(i)))), trim(dist_names(dist_index(i))), &
                  trim(lambda_names(lambda_index(i))), results(i)%lambda_on, results(i)%theta_on, results(i)%gamma_on, &
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
        print '(A18,1X,A8,1X,A10,1X,A10,1X,A10,1X,A4,1X,A12,1X,A12,1X,A12,1X,A12)', &
              "Model", "Dist", "lambda_on", "theta_on", "gamma_on", "df", "LR_vs_ON0", "p_value", &
              "dAIC_ONEST", "dBIC_ONEST"
        print '(A)', "----------------------------------------------------------------------------------------------------"
        do idist = 1, size(dist_names)
            do imodel = 1, size(model_names)
                idx0 = flat_index(imodel, idist, lambda_zero)
                idxest = flat_index(imodel, idist, lambda_est)
                df = lambda_param_count(model_id(trim(model_names(imodel))), lambda_est)
                lr_stat = max(2.0_dp * (results(idxest)%fit%loglik - results(idx0)%fit%loglik), 0.0_dp)
                if (df == 2) then
                    p_value = exp(-0.5_dp * lr_stat)
                else
                    p_value = erfc(sqrt(0.5_dp * lr_stat))
                end if
                delta_aic = results(idxest)%aic - results(idx0)%aic
                delta_bic = results(idxest)%bic - results(idx0)%bic
                print '(A18,1X,A8,1X,F10.4,1X,F10.4,1X,F10.4,1X,I4,1X,F12.3,1X,ES12.4,1X,F12.3,1X,F12.3)', &
                      model_table_label(trim(model_names(imodel))), trim(dist_names(idist)), &
                      results(idxest)%lambda_on, results(idxest)%theta_on, results(idxest)%gamma_on, df, &
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

    ! Dispatch to the selected public MCS-GARCH filter.
    subroutine filter_by_model(model, y, daily_var, diurnal_var, params, q)
        integer, intent(in) :: model
        real(dp), intent(in) :: y(:), daily_var(:), diurnal_var(:)
        type(mcsgarch_params_t), intent(in) :: params
        real(dp), intent(out) :: q(:)

        select case (model)
        case (model_nagarch)
            call mcsgarch_nagarch_filter(y, daily_var, diurnal_var, params, q)
        case (model_gjr)
            call mcsgarch_gjr_filter(y, daily_var, diurnal_var, params, q)
        case default
            call mcsgarch_filter(y, daily_var, diurnal_var, params, q)
        end select
    end subroutine filter_by_model

    ! Build the intraday daily scale using an overnight NIC matching the selected intraday model family.
    subroutine intraday_scale_from_overnight(base_scale, overnight_return, model, lambda_on, theta_on, gamma_on, scale)
        real(dp), intent(in) :: base_scale(:), overnight_return(:), lambda_on, theta_on, gamma_on
        integer, intent(in) :: model
        real(dp), intent(out) :: scale(:)
        integer :: i
        real(dp) :: r, impact

        if (size(base_scale) /= size(overnight_return) .or. size(scale) /= size(base_scale)) then
            error stop "intraday_scale_from_overnight: array sizes differ"
        end if
        do i = 1, size(base_scale)
            r = overnight_return(i)
            select case (model)
            case (model_nagarch)
                impact = lambda_on * (r - theta_on * sqrt(max(base_scale(i), min_var)))**2
            case (model_gjr)
                impact = lambda_on * r**2 + gamma_on * merge(r**2, 0.0_dp, r < 0.0_dp)
            case default
                impact = lambda_on * r**2
            end select
            scale(i) = max(base_scale(i) + impact, min_var)
        end do
    end subroutine intraday_scale_from_overnight

    ! Pack constrained model, overnight, and distribution parameters.
    subroutine pack_joint_params(params, overnight_weight, lambda_on, theta_on, gamma_on, model, dist, lambda_mode, p)
        type(mcsgarch_params_t), intent(in) :: params
        real(dp), intent(in) :: overnight_weight, lambda_on, theta_on, gamma_on
        integer, intent(in) :: model, dist, lambda_mode
        real(dp), intent(out) :: p(:)
        real(dp) :: rest, a, b, ghalf, theta
        integer :: j

        p(1) = log(max(params%omega, min_var))
        theta = params%theta
        if (model == model_nagarch) then
            a = max(params%alpha * (1.0_dp + theta**2), min_var)
        else
            a = max(params%alpha, min_var)
        end if
        b = max(params%beta, min_var)
        if (model == model_gjr) then
            ghalf = max(0.5_dp*params%gamma, min_var)
            rest = max(persist_max - a - b - ghalf, min_var)
        else
            ghalf = 0.0_dp
            rest = max(persist_max - a - b, min_var)
        end if
        p(2) = log(a / rest)
        p(3) = log(b / rest)
        if (model == model_nagarch) p(4) = theta
        if (model == model_gjr) p(4) = log(ghalf / rest)
        j = model_param_count(model) + 1
        p(j) = log(max(overnight_weight, min_var))
        if (lambda_mode == lambda_est) then
            j = j + 1
            p(j) = log(max(lambda_on, min_var))
            if (model == model_nagarch) then
                j = j + 1
                p(j) = theta_on
            else if (model == model_gjr) then
                j = j + 1
                p(j) = log(max(gamma_on, min_var))
            end if
        end if
        if (dist == dist_t) p(j + 1) = log((8.0_dp - 2.0_dp) / (100.0_dp - 8.0_dp))
        if (dist == dist_fs_skewt) then
            p(j + 1) = log((8.0_dp - 2.0_dp) / (100.0_dp - 8.0_dp))
            p(j + 2) = 0.0_dp
        end if
    end subroutine pack_joint_params

    ! Unpack optimizer variables into constrained parameters.
    subroutine unpack_joint_params(p, model, dist, params, overnight_weight, lambda_on, theta_on, gamma_on, nu, xi, &
                                   lambda_mode)
        real(dp), intent(in) :: p(:)
        integer, intent(in) :: model, dist, lambda_mode
        type(mcsgarch_params_t), intent(out) :: params
        real(dp), intent(out) :: overnight_weight, lambda_on, theta_on, gamma_on, nu, xi
        real(dp) :: ea, eb, eg, den, theta, alpha_effective
        integer :: j

        params%omega = exp(max(min(p(1), 50.0_dp), -50.0_dp))
        ea = exp(max(min(p(2), 50.0_dp), -50.0_dp))
        eb = exp(max(min(p(3), 50.0_dp), -50.0_dp))
        if (model == model_gjr) then
            eg = exp(max(min(p(4), 50.0_dp), -50.0_dp))
        else
            eg = 0.0_dp
        end if
        den = 1.0_dp + ea + eb + eg
        alpha_effective = persist_max * ea / den
        params%beta = persist_max * eb / den
        params%gamma = 0.0_dp
        if (model == model_nagarch) then
            theta = max(min(p(4), 20.0_dp), -20.0_dp)
            params%theta = theta
            params%alpha = alpha_effective / (1.0_dp + theta**2)
        else if (model == model_gjr) then
            params%theta = 0.0_dp
            params%alpha = alpha_effective
            params%gamma = 2.0_dp * persist_max * eg / den
        else
            params%theta = 0.0_dp
            params%alpha = alpha_effective
        end if
        j = model_param_count(model) + 1
        overnight_weight = exp(max(min(p(j), 50.0_dp), -50.0_dp))
        select case (lambda_mode)
        case (lambda_zero)
            lambda_on = 0.0_dp
            theta_on = 0.0_dp
            gamma_on = 0.0_dp
        case (lambda_est)
            j = j + 1
            lambda_on = exp(max(min(p(j), 50.0_dp), -50.0_dp))
            theta_on = 0.0_dp
            gamma_on = 0.0_dp
            if (model == model_nagarch) then
                j = j + 1
                theta_on = max(min(p(j), 20.0_dp), -20.0_dp)
            else if (model == model_gjr) then
                j = j + 1
                gamma_on = exp(max(min(p(j), 50.0_dp), -50.0_dp))
            end if
        case (lambda_one)
            lambda_on = 1.0_dp
            theta_on = 0.0_dp
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

    ! Build initial dynamic parameters from scaled intraday residual variance.
    type(mcsgarch_params_t) function initial_params(intra_ret, intra_scale, diurnal_var)
        real(dp), intent(in) :: intra_ret(:), intra_scale(:), diurnal_var(:)
        real(dp) :: vnorm

        vnorm = max(variance(intra_ret / sqrt(max(intra_scale * diurnal_var, min_var))), min_var)
        initial_params = mcsgarch_params_t(omega=vnorm * 0.04_dp, alpha=0.06_dp, beta=0.90_dp, &
                                           theta=0.0_dp, gamma=0.05_dp)
    end function initial_params

    ! Estimate the starting overnight variance weight by method of moments.
    real(dp) function initial_overnight_weight(on_ret, on_scale)
        real(dp), intent(in) :: on_ret(:), on_scale(:)

        initial_overnight_weight = max(mean(on_ret**2 / max(on_scale, min_var)), min_var)
    end function initial_overnight_weight

    ! Return whether parameters satisfy positivity and persistence constraints.
    logical function params_valid(params, model)
        type(mcsgarch_params_t), intent(in) :: params
        integer, intent(in) :: model

        params_valid = params%omega > 0.0_dp .and. params%alpha >= 0.0_dp .and. &
                       params%beta >= 0.0_dp .and. params%gamma >= 0.0_dp .and. &
                       persist(params, model) < 1.0_dp
    end function params_valid

    ! Return persistence measure for the selected q_t dynamics.
    real(dp) function persist(params, model)
        type(mcsgarch_params_t), intent(in) :: params
        integer, intent(in) :: model

        select case (model)
        case (model_nagarch)
            persist = params%alpha * (1.0_dp + params%theta**2) + params%beta
        case (model_gjr)
            persist = params%alpha + 0.5_dp*params%gamma + params%beta
        case default
            persist = params%alpha + params%beta
        end select
    end function persist

    ! Return model id for a configured model name.
    integer function model_id(model_name)
        character(len=*), intent(in) :: model_name

        select case (trim(model_name))
        case ("MCSGARCH")
            model_id = model_sym
        case ("MCSNAGARCH")
            model_id = model_nagarch
        case ("MCSGJRGARCH")
            model_id = model_gjr
        case default
            error stop "model_id: unsupported model"
        end select
    end function model_id

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

    ! Return the dynamic-parameter count for an MCS model id.
    integer function model_param_count(model)
        integer, intent(in) :: model

        if (model == model_sym) then
            model_param_count = 3
        else
            model_param_count = 4
        end if
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
    integer function lambda_param_count(model, lambda_mode)
        integer, intent(in) :: model, lambda_mode

        if (lambda_mode == lambda_est) then
            if (model == model_sym) then
                lambda_param_count = 1
            else
                lambda_param_count = 2
            end if
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
    real(dp) function initial_gamma_on(model, lambda_mode)
        integer, intent(in) :: model, lambda_mode

        if (model == model_gjr .and. lambda_mode == lambda_est) then
            initial_gamma_on = 0.05_dp
        else
            initial_gamma_on = 0.0_dp
        end if
    end function initial_gamma_on

    ! Return display label used in output tables.
    character(len=36) function model_table_label(model_name)
        character(len=*), intent(in) :: model_name

        select case (trim(model_name))
        case ("MCSGARCH")
            model_table_label = "MCS-GARCH"
        case ("MCSNAGARCH")
            model_table_label = "MCS-NAGARCH"
        case ("MCSGJRGARCH")
            model_table_label = "MCS-GJRGARCH"
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

end module fit_mcsgarch_on_intraday_mod

program xfit_mcsgarch_on_intraday
    use date_mod, only: print_program_header
    use fit_mcsgarch_on_intraday_mod, only: run_fit_mcsgarch_on_intraday
    implicit none
    call print_program_header("xfit_mcsgarch_on_intraday.f90")
    call run_fit_mcsgarch_on_intraday()
end program xfit_mcsgarch_on_intraday

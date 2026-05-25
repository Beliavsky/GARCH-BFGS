program xsim_garch_fit
    use, intrinsic :: iso_fortran_env, only: int64
    use kind_mod, only: dp
    use garch_types_mod, only: garch_params_t
    use garch_sim_mod, only: simulate_garch_model
    use garch_fit_mod, only: fit_symm_garch, fit_nagarch, fit_gjr, fit_gjr_signed, fit_egarch, &
                             fit_qgarch, fit_csgarch, fit_tgarch, &
                             fit_symm_garch_pq, fit_nagarch_pq, fit_figarch, fit_fi_nagarch, &
                             fit_harch, fit_avgarch, fit_aparch, fit_midas_hyperbolic, &
                             fit_midas_hyperbolic_asym, fit_fgarch_twist, fit_ewma, &
                             symm_garch_persist, nagarch_persist, gjr_persist, egarch_persist, &
                             qgarch_persist, csgarch_persist, tgarch_persist, &
                             symm_garch_pq_persist, nagarch_pq_persist, figarch_persist, &
                             fi_nagarch_persist, harch_persist, avgarch_persist, aparch_persist, &
                             midas_hyperbolic_persist, midas_hyperbolic_asym_persist, fgarch_twist_persist
    implicit none

    integer, parameter :: n_model = 18
    integer, parameter :: n_keep = 3000
    integer, parameter :: n_burn = 500
    integer, parameter :: max_iter = 400
    integer, parameter :: n_stat = 6
    integer, parameter :: n_oos = 4
    real(dp), parameter :: gtol = 1.0e-6_dp
    real(dp), parameter :: rm_lambda = 0.94_dp

    character(len=16), parameter :: model_names(n_model) = [ &
        "SYMM_GARCH      ", &
        "SYMM_GARCH21    ", &
        "NAGARCH         ", &
        "NAGARCH21       ", &
        "GJR_GARCH       ", &
        "GJR_SIGNED      ", &
        "EGARCH          ", &
        "QGARCH          ", &
        "CSGARCH         ", &
        "TGARCH          ", &
        "FIGARCH         ", &
        "FI_NAGARCH      ", &
        "HARCH           ", &
        "AVGARCH         ", &
        "APARCH          ", &
        "MIDAS_HYPER     ", &
        "MIDAS_ASYM      ", &
        "FGARCH_TWIST    " ]
    character(len=7), parameter :: stat_names(n_stat) = [ &
        "omega  ", "alpha  ", "gamma  ", "beta   ", "theta  ", "persist" ]
    character(len=12), parameter :: oos_names(n_oos) = [ &
        "TRUE_FIT    ", "RM_EWMA     ", "FIT_EWMA    ", "GARCH11     " ]

    type(garch_params_t) :: true_params, fit_params
    real(dp), allocatable :: y_all(:), y(:), y_train_test(:)
    real(dp) :: f_best
    real(dp) :: true_persist, fit_persist
    real(dp) :: true_vals(n_stat), fit_vals(n_stat)
    real(dp) :: sum_fit(n_stat), sum_sqerr(n_stat)
    real(dp) :: oos_loss(n_oos), sum_oos(n_oos), win_oos(n_oos)
    real(dp) :: oos_win_rate_model(n_model,n_oos)
    real(dp) :: sim_time_model(n_model), fit_time_model(n_model)
    real(dp) :: sum_sim_time, sum_fit_time
    integer :: i, irep, niter, seed_val, n_rep, ntest, n_total, ios
    integer(int64) :: clock_count0, clock_count1, clock_rate
    integer :: conv_count
    character(len=64) :: arg1
    logical :: converged

    n_rep = 1
    ntest = 0
    call system_clock(count_rate=clock_rate)
    oos_win_rate_model = 0.0_dp
    sim_time_model = 0.0_dp
    fit_time_model = 0.0_dp
    if (command_argument_count() >= 1) then
        call get_command_argument(1, arg1)
        read(arg1, *, iostat=ios) n_rep
        if (ios /= 0 .or. n_rep < 1) then
            write(*,'(A)') "Usage: xsim_garch_fit.exe [n_rep>=1] [ntest>=0]"
            error stop "xsim_garch_fit: invalid n_rep argument"
        end if
    end if
    if (command_argument_count() >= 2) then
        call get_command_argument(2, arg1)
        read(arg1, *, iostat=ios) ntest
        if (ios /= 0 .or. ntest < 0) then
            write(*,'(A)') "Usage: xsim_garch_fit.exe [n_rep>=1] [ntest>=0]"
            error stop "xsim_garch_fit: invalid ntest argument"
        end if
    end if

    n_total = n_keep + n_burn + ntest
    allocate(y_all(n_total), y(n_keep))
    if (ntest > 0) then
        allocate(y_train_test(n_keep + ntest))
    end if

    write(*,'(A)') "Simulate -> fit sanity check"
    write(*,'(A,I0)') "Observations simulated per dataset (incl burn-in): ", n_total
    write(*,'(A,I0,A,I0,A,I0,A,I0)') "Kept observations: ", n_keep, &
        "  burn-in: ", n_burn, "  datasets/model: ", n_rep, "  models: ", n_model
    if (ntest > 0) then
        write(*,'(A,I0,A)') "Out-of-sample observations: ", ntest, " (enabled)"
    else
        write(*,'(A,I0,A)') "Out-of-sample observations: ", ntest, " (disabled)"
    end if
    write(*,'(A)') "Rows: TRUE parameters and FIT estimates from simulated returns"
    write(*,'(A)') repeat("-", 146)
    write(*,'(A)') "Model            Row   Set      omega    alpha    gamma     beta    theta  persist  niter conv    nll/obs"
    write(*,'(A)') repeat("-", 146)

    do i = 1, n_model
        call set_true_params(trim(model_names(i)), true_params)
        true_persist = model_persist(trim(model_names(i)), true_params)
        call collect_stats(true_params, true_persist, true_vals)
        sum_fit = 0.0_dp
        sum_sqerr = 0.0_dp
        sum_oos = 0.0_dp
        win_oos = 0.0_dp
        sum_sim_time = 0.0_dp
        sum_fit_time = 0.0_dp
        conv_count = 0

        write(*,'(A16,1X,A5,1X,A4,1X,ES9.2,1X,F8.4,1X,F8.4,1X,F8.4,1X,F8.4,1X,F8.4,1X,A5,1X,A4,1X,A10)') &
            trim(model_names(i)), "TRUE", "-", true_params%omega, true_params%alpha, true_params%gamma, &
            true_params%beta, true_params%theta, true_persist, "-", "-", "-"

        do irep = 1, n_rep
            seed_val = 100000 + 1000*i + irep
            call system_clock(clock_count0)
            call simulate_garch_model(trim(model_names(i)), true_params, n_total, seed_val, y_all)
            call system_clock(clock_count1)
            sum_sim_time = sum_sim_time + elapsed_seconds(clock_count0, clock_count1, clock_rate)
            y = y_all(n_burn + 1:n_burn + n_keep)
            if (ntest > 0) then
                y_train_test = y_all(n_burn + 1:n_total)
            end if

            call system_clock(clock_count0)
            call fit_one_model(trim(model_names(i)), y, max_iter, gtol, f_best, fit_params, niter, converged)
            call system_clock(clock_count1)
            sum_fit_time = sum_fit_time + elapsed_seconds(clock_count0, clock_count1, clock_rate)
            fit_persist = model_persist(trim(model_names(i)), fit_params)
            call collect_stats(fit_params, fit_persist, fit_vals)
            sum_fit = sum_fit + fit_vals
            sum_sqerr = sum_sqerr + (fit_vals - true_vals)**2
            if (converged) conv_count = conv_count + 1
            if (ntest > 0) then
                call score_oos(trim(model_names(i)), y_train_test, n_keep, ntest, max_iter, gtol, &
                               fit_params, oos_loss)
                sum_oos = sum_oos + oos_loss
                where (oos_loss(1) < oos_loss(2:n_oos))
                    win_oos(2:n_oos) = win_oos(2:n_oos) + 1.0_dp
                end where
            end if

            write(*,'(A16,1X,A5,1X,I4,1X,ES9.2,1X,F8.4,1X,F8.4,1X,F8.4,1X,F8.4,1X,F8.4,1X,I5,1X,L4,1X,F10.6)') &
                trim(model_names(i)), "FIT", irep, fit_params%omega, fit_params%alpha, fit_params%gamma, &
                fit_params%beta, fit_params%theta, fit_persist, niter, converged, f_best
        end do
        call print_summary(trim(model_names(i)), true_vals, sum_fit, sum_sqerr, n_rep, conv_count)
        sim_time_model(i) = sum_sim_time / real(n_rep, dp)
        fit_time_model(i) = sum_fit_time / real(n_rep, dp)
        if (ntest > 0) then
            oos_win_rate_model(i,:) = win_oos / real(n_rep, dp)
            call print_oos_summary(trim(model_names(i)), sum_oos, win_oos, n_rep)
        end if
        write(*,'(A)') repeat("-", 146)
    end do

    if (ntest > 0) call print_oos_win_table(oos_win_rate_model, sim_time_model, fit_time_model)

contains

    real(dp) function elapsed_seconds(count0, count1, count_rate)
        integer(int64), intent(in) :: count0, count1, count_rate

        if (count_rate <= 0_int64 .or. count1 < count0) then
            elapsed_seconds = 0.0_dp
        else
            elapsed_seconds = real(count1 - count0, dp) / real(count_rate, dp)
        end if
    end function elapsed_seconds

    subroutine collect_stats(params, persist, vals)
        type(garch_params_t), intent(in) :: params
        real(dp), intent(in) :: persist
        real(dp), intent(out) :: vals(n_stat)

        vals = [params%omega, params%alpha, params%gamma, params%beta, params%theta, persist]
    end subroutine collect_stats

    subroutine print_summary(model_name, true_vals, sum_fit, sum_sqerr, n_rep, conv_count)
        character(len=*), intent(in) :: model_name
        real(dp), intent(in) :: true_vals(n_stat), sum_fit(n_stat), sum_sqerr(n_stat)
        integer, intent(in) :: n_rep, conv_count
        real(dp) :: mean_fit(n_stat), bias(n_stat), rmse(n_stat), conv_rate
        integer :: j

        mean_fit = sum_fit / real(n_rep, dp)
        bias = mean_fit - true_vals
        rmse = sqrt(sum_sqerr / real(n_rep, dp))
        conv_rate = real(conv_count, dp) / real(n_rep, dp)

        write(*,'(A,1X,A,1X,A,I0,A,I0,A,F7.3)') &
            "SUMMARY", trim(model_name), "converged=", conv_count, "/", n_rep, " rate=", conv_rate
        write(*,'(A)') "          Param          true      mean_fit          bias          rmse"
        do j = 1, n_stat
            write(*,'(10X,A7,1X,ES12.4,1X,ES12.4,1X,ES12.4,1X,ES12.4)') &
                stat_names(j), true_vals(j), mean_fit(j), bias(j), rmse(j)
        end do
    end subroutine print_summary

    subroutine print_oos_summary(model_name, sum_oos, win_oos, n_rep)
        character(len=*), intent(in) :: model_name
        real(dp), intent(in) :: sum_oos(n_oos), win_oos(n_oos)
        integer, intent(in) :: n_rep
        real(dp) :: mean_oos(n_oos)
        integer :: j

        mean_oos = sum_oos / real(n_rep, dp)
        write(*,'(A,1X,A)') "OOS_QLIKE", trim(model_name)
        write(*,'(A)') "          Model        mean_qlike  true_fit_win_rate"
        write(*,'(10X,A12,1X,F11.6,1X,A)') oos_names(1), mean_oos(1), "-"
        do j = 2, n_oos
            write(*,'(10X,A12,1X,F11.6,1X,F17.3)') &
                oos_names(j), mean_oos(j), win_oos(j) / real(n_rep, dp)
        end do
    end subroutine print_oos_summary

    subroutine print_oos_win_table(oos_win_rate_model, sim_time_model, fit_time_model)
        real(dp), intent(in) :: oos_win_rate_model(n_model,n_oos)
        real(dp), intent(in) :: sim_time_model(n_model), fit_time_model(n_model)
        integer :: im

        write(*,'(/,A)') "OOS true_fit_win_rate by simulated model"
        write(*,'(A)') repeat("-", 108)
        write(*,'(A16,1X,A11,1X,A11,1X,A11,1X,A12,1X,A12)') &
            "Model", "vs_RM_EWMA", "vs_FIT_EWMA", "vs_GARCH11", "sim_sec", "fit_sec"
        write(*,'(A)') repeat("-", 108)
        do im = 1, n_model
            write(*,'(A16,1X,F11.3,1X,F11.3,1X,F11.3,1X,F12.6,1X,F12.6)') trim(model_names(im)), &
                oos_win_rate_model(im,2), oos_win_rate_model(im,3), oos_win_rate_model(im,4), &
                sim_time_model(im), fit_time_model(im)
        end do
        write(*,'(A)') repeat("-", 108)
    end subroutine print_oos_win_table

    subroutine score_oos(model_name, y_full, ntrain, ntest, max_iter, gtol, true_fit_params, loss)
        character(len=*), intent(in) :: model_name
        real(dp), intent(in) :: y_full(:), gtol
        integer, intent(in) :: ntrain, ntest, max_iter
        type(garch_params_t), intent(in) :: true_fit_params
        real(dp), intent(out) :: loss(n_oos)
        type(garch_params_t) :: ewma_params, garch_params
        real(dp) :: f_best
        integer :: niter
        logical :: converged

        loss(1) = oos_qlike_model(model_name, y_full, ntrain, ntest, true_fit_params)

        ewma_params = garch_params_t()
        ewma_params%beta = rm_lambda
        ewma_params%alpha = 1.0_dp - rm_lambda
        loss(2) = oos_qlike_model("EWMA", y_full, ntrain, ntest, ewma_params)

        call fit_ewma(y_full(1:ntrain), max_iter, gtol, f_best, ewma_params, niter, converged)
        loss(3) = oos_qlike_model("EWMA", y_full, ntrain, ntest, ewma_params)

        call fit_symm_garch(y_full(1:ntrain), max_iter, gtol, f_best, garch_params, niter, converged)
        loss(4) = oos_qlike_model("SYMM_GARCH", y_full, ntrain, ntest, garch_params)
    end subroutine score_oos

    real(dp) function oos_qlike_model(model_name, y_full, ntrain, ntest, params)
        character(len=*), intent(in) :: model_name
        real(dp), intent(in) :: y_full(:)
        integer, intent(in) :: ntrain, ntest
        type(garch_params_t), intent(in) :: params
        real(dp), allocatable :: h(:)
        integer :: t

        allocate(h(size(y_full)))
        call model_variance_path(model_name, y_full, ntrain, params, h)
        oos_qlike_model = 0.0_dp
        do t = ntrain + 1, ntrain + ntest
            h(t) = max(h(t), 1.0e-12_dp)
            oos_qlike_model = oos_qlike_model + log(h(t)) + y_full(t)**2 / h(t)
        end do
        oos_qlike_model = oos_qlike_model / real(ntest, dp)
        deallocate(h)
    end function oos_qlike_model

    subroutine model_variance_path(model_name, y_full, ntrain, params, h)
        character(len=*), intent(in) :: model_name
        real(dp), intent(in) :: y_full(:)
        integer, intent(in) :: ntrain
        type(garch_params_t), intent(in) :: params
        real(dp), intent(out) :: h(:)

        select case (trim(adjustl(model_name)))
        case ("EWMA")
            call ewma_variance_path(y_full, ntrain, params%beta, h)
        case ("SYMM_GARCH", "GARCH")
            call garch11_variance_path(y_full, ntrain, params, h)
        case ("NAGARCH")
            call nagarch_variance_path(y_full, ntrain, params, h)
        case ("GJR_GARCH", "GJR", "GJR_SIGNED")
            call gjr_variance_path(y_full, ntrain, params, h)
        case ("EGARCH")
            call egarch_variance_path(y_full, params, h)
        case ("QGARCH")
            call qgarch_variance_path(y_full, params, h)
        case ("CSGARCH")
            call csgarch_variance_path(y_full, ntrain, params, h)
        case ("TGARCH")
            call tgarch_variance_path(y_full, ntrain, params, h)
        case ("SYMM_GARCH21", "SYMM_GARCH_PQ")
            call garch_pq_variance_path(y_full, ntrain, params, h)
        case ("NAGARCH21", "NAGARCH_PQ")
            call nagarch_pq_variance_path(y_full, ntrain, params, h)
        case ("FIGARCH")
            call figarch_variance_path(y_full, ntrain, params, h)
        case ("FI_NAGARCH")
            call fi_nagarch_variance_path(y_full, ntrain, params, h)
        case ("HARCH")
            call harch_variance_path(y_full, ntrain, params, h)
        case ("AVGARCH")
            call avgarch_variance_path(y_full, ntrain, params, h)
        case ("APARCH")
            call aparch_variance_path(y_full, ntrain, params, h)
        case ("MIDAS_HYPER", "MIDAS_HYPERBOLIC")
            call midas_variance_path(y_full, ntrain, params, h)
        case ("MIDAS_ASYM", "MIDAS_HYPERBOLIC_ASYM")
            call midas_asym_variance_path(y_full, ntrain, params, h)
        case ("FGARCH_TWIST")
            call fgarch_twist_variance_path(y_full, params, h)
        case default
            error stop "xsim_garch_fit: no OOS variance path for model"
        end select
    end subroutine model_variance_path

    subroutine ewma_variance_path(y, ntrain, lambda, h)
        real(dp), intent(in) :: y(:), lambda
        integer, intent(in) :: ntrain
        real(dp), intent(out) :: h(:)
        real(dp) :: ht
        integer :: t

        ht = sample_variance_n(y, ntrain)
        do t = 1, size(y)
            h(t) = max(ht, 1.0e-12_dp)
            ht = lambda*h(t) + (1.0_dp - lambda)*y(t)**2
        end do
    end subroutine ewma_variance_path

    subroutine garch11_variance_path(y, ntrain, params, h)
        real(dp), intent(in) :: y(:)
        integer, intent(in) :: ntrain
        type(garch_params_t), intent(in) :: params
        real(dp), intent(out) :: h(:)
        integer :: t

        h(1) = max(params%omega / max(1.0_dp - params%alpha - params%beta, 1.0e-8_dp), &
                   sample_variance_n(y, ntrain))
        h(1) = max(h(1), 1.0e-12_dp)
        do t = 2, size(y)
            h(t) = params%omega + params%alpha*y(t - 1)**2 + params%beta*h(t - 1)
            h(t) = max(h(t), 1.0e-12_dp)
        end do
    end subroutine garch11_variance_path

    subroutine nagarch_variance_path(y, ntrain, params, h)
        real(dp), intent(in) :: y(:)
        integer, intent(in) :: ntrain
        type(garch_params_t), intent(in) :: params
        real(dp), intent(out) :: h(:)
        real(dp) :: sqrth
        integer :: t

        h(1) = max(params%omega / max(1.0_dp - nagarch_persist(params), 1.0e-8_dp), &
                   sample_variance_n(y, ntrain))
        h(1) = max(h(1), 1.0e-12_dp)
        do t = 2, size(y)
            sqrth = sqrt(max(h(t - 1), 1.0e-12_dp))
            h(t) = params%omega + params%alpha*(y(t - 1) - params%theta*sqrth)**2 + params%beta*h(t - 1)
            h(t) = max(h(t), 1.0e-12_dp)
        end do
    end subroutine nagarch_variance_path

    subroutine gjr_variance_path(y, ntrain, params, h)
        real(dp), intent(in) :: y(:)
        integer, intent(in) :: ntrain
        type(garch_params_t), intent(in) :: params
        real(dp), intent(out) :: h(:)
        real(dp) :: ind
        integer :: t

        h(1) = max(params%omega / max(1.0_dp - gjr_persist(params), 1.0e-8_dp), &
                   sample_variance_n(y, ntrain))
        h(1) = max(h(1), 1.0e-12_dp)
        do t = 2, size(y)
            ind = merge(1.0_dp, 0.0_dp, y(t - 1) < 0.0_dp)
            h(t) = params%omega + (params%alpha + params%gamma*ind)*y(t - 1)**2 + params%beta*h(t - 1)
            h(t) = max(h(t), 1.0e-12_dp)
        end do
    end subroutine gjr_variance_path

    subroutine egarch_variance_path(y, params, h)
        real(dp), intent(in) :: y(:)
        type(garch_params_t), intent(in) :: params
        real(dp), intent(out) :: h(:)
        real(dp) :: lh, z, c_eg
        integer :: t

        c_eg = sqrt(2.0_dp / acos(-1.0_dp))
        lh = params%omega / max(1.0_dp - params%beta, 1.0e-8_dp)
        do t = 1, size(y)
            h(t) = max(exp(lh), 1.0e-12_dp)
            z = y(t) / sqrt(h(t))
            lh = params%omega + params%beta*lh + params%alpha*(abs(z) - c_eg) + params%gamma*z
        end do
    end subroutine egarch_variance_path

    subroutine qgarch_variance_path(y, params, h)
        real(dp), intent(in) :: y(:)
        type(garch_params_t), intent(in) :: params
        real(dp), intent(out) :: h(:)
        integer :: t

        h(1) = max((params%omega + params%alpha*params%theta**2) / &
                   max(1.0_dp - params%alpha - params%beta, 1.0e-8_dp), 1.0e-12_dp)
        do t = 2, size(y)
            h(t) = params%omega + params%alpha*(y(t - 1) - params%theta)**2 + params%beta*h(t - 1)
            h(t) = max(h(t), 1.0e-12_dp)
        end do
    end subroutine qgarch_variance_path

    subroutine csgarch_variance_path(y, ntrain, params, h)
        real(dp), intent(in) :: y(:)
        integer, intent(in) :: ntrain
        type(garch_params_t), intent(in) :: params
        real(dp), intent(out) :: h(:)
        real(dp) :: q_prev, h_prev, shock, qq, backcast
        integer :: t

        backcast = sample_variance_n(y, ntrain)
        q_prev = max(params%omega / max(1.0_dp - params%extra1, 1.0e-8_dp), 1.0e-12_dp)
        h_prev = max(q_prev, backcast)
        shock = backcast
        do t = 1, size(y)
            qq = max(params%omega + params%extra1*q_prev + params%extra2*(shock - h_prev), 1.0e-12_dp)
            h(t) = max(qq + params%alpha*(shock - q_prev) + params%beta*(h_prev - q_prev), 1.0e-12_dp)
            q_prev = qq
            h_prev = h(t)
            shock = y(t)**2
        end do
    end subroutine csgarch_variance_path

    subroutine tgarch_variance_path(y, ntrain, params, h)
        real(dp), intent(in) :: y(:)
        integer, intent(in) :: ntrain
        type(garch_params_t), intent(in) :: params
        real(dp), intent(out) :: h(:)
        real(dp) :: sigma, z, news
        integer :: t

        sigma = max(params%omega / max(1.0_dp - tgarch_persist(params), 1.0e-8_dp), &
                    sqrt(sample_variance_n(y, ntrain)))
        do t = 1, size(y)
            sigma = max(sigma, 1.0e-6_dp)
            h(t) = sigma**2
            z = y(t) / sigma
            news = sqrt(1.0e-6_dp + z**2) - params%gamma*z
            sigma = params%omega + params%alpha*news*sigma + params%beta*sigma
        end do
    end subroutine tgarch_variance_path

    subroutine garch_pq_variance_path(y, ntrain, params, h)
        real(dp), intent(in) :: y(:)
        integer, intent(in) :: ntrain
        type(garch_params_t), intent(in) :: params
        real(dp), intent(out) :: h(:)
        real(dp) :: backcast, ht, lag_h
        integer :: t, i, lag

        backcast = sample_variance_n(y, ntrain)
        do t = 1, size(y)
            ht = params%omega
            do i = 1, size(params%alpha_lags)
                lag = t - i
                if (lag >= 1) then
                    ht = ht + params%alpha_lags(i)*y(lag)**2
                else
                    ht = ht + params%alpha_lags(i)*backcast
                end if
            end do
            do i = 1, size(params%beta_lags)
                lag = t - i
                if (lag >= 1) then
                    lag_h = h(lag)
                else
                    lag_h = backcast
                end if
                ht = ht + params%beta_lags(i)*lag_h
            end do
            h(t) = max(ht, 1.0e-12_dp)
        end do
    end subroutine garch_pq_variance_path

    subroutine nagarch_pq_variance_path(y, ntrain, params, h)
        real(dp), intent(in) :: y(:)
        integer, intent(in) :: ntrain
        type(garch_params_t), intent(in) :: params
        real(dp), intent(out) :: h(:)
        real(dp) :: backcast, ht, lag_h, lag_y, sqrth
        integer :: t, i, lag

        backcast = sample_variance_n(y, ntrain)
        do t = 1, size(y)
            ht = params%omega
            do i = 1, size(params%alpha_lags)
                lag = t - i
                if (lag >= 1) then
                    lag_h = h(lag)
                    lag_y = y(lag)
                else
                    lag_h = backcast
                    lag_y = 0.0_dp
                end if
                sqrth = sqrt(max(lag_h, 1.0e-12_dp))
                ht = ht + params%alpha_lags(i)*(lag_y - params%theta*sqrth)**2
            end do
            do i = 1, size(params%beta_lags)
                lag = t - i
                if (lag >= 1) then
                    lag_h = h(lag)
                else
                    lag_h = backcast
                end if
                ht = ht + params%beta_lags(i)*lag_h
            end do
            h(t) = max(ht, 1.0e-12_dp)
        end do
    end subroutine nagarch_pq_variance_path

    subroutine figarch_variance_path(y, ntrain, params, h)
        real(dp), intent(in) :: y(:)
        integer, intent(in) :: ntrain
        type(garch_params_t), intent(in) :: params
        real(dp), intent(out) :: h(:)
        real(dp), allocatable :: lambda(:)
        real(dp) :: omega_tilde, backcast, bc_weight
        integer :: t, i, m

        m = 1000
        allocate(lambda(m))
        call figarch_weights_local(params%alpha, params%theta, params%beta, lambda)
        omega_tilde = params%omega / max(1.0_dp - params%beta, 1.0e-8_dp)
        backcast = sample_variance_n(y, ntrain)
        do t = 1, size(y)
            bc_weight = 0.0_dp
            do i = t, m
                bc_weight = bc_weight + lambda(i)
            end do
            h(t) = omega_tilde + bc_weight*backcast
            do i = 1, min(t - 1, m)
                h(t) = h(t) + lambda(i)*y(t - i)**2
            end do
            h(t) = max(h(t), 1.0e-12_dp)
        end do
        deallocate(lambda)
    end subroutine figarch_variance_path

    subroutine fi_nagarch_variance_path(y, ntrain, params, h)
        real(dp), intent(in) :: y(:)
        integer, intent(in) :: ntrain
        type(garch_params_t), intent(in) :: params
        real(dp), intent(out) :: h(:)
        real(dp), allocatable :: lambda(:), news(:)
        real(dp) :: omega_tilde, backcast, bc_weight, sqrth, scale
        integer :: t, i, m

        m = 1000
        allocate(lambda(m), news(size(y)))
        call figarch_weights_local(params%alpha, params%theta, params%beta, lambda)
        omega_tilde = params%omega / max(1.0_dp - params%beta, 1.0e-8_dp)
        backcast = sample_variance_n(y, ntrain)
        scale = 1.0_dp + params%twist**2
        do t = 1, size(y)
            bc_weight = 0.0_dp
            do i = t, m
                bc_weight = bc_weight + lambda(i)
            end do
            h(t) = omega_tilde + bc_weight*backcast
            do i = 1, min(t - 1, m)
                h(t) = h(t) + lambda(i)*news(t - i)
            end do
            h(t) = max(h(t), 1.0e-12_dp)
            sqrth = sqrt(h(t))
            news(t) = (y(t) - params%twist*sqrth)**2 / scale
        end do
        deallocate(lambda, news)
    end subroutine fi_nagarch_variance_path

    subroutine harch_variance_path(y, ntrain, params, h)
        real(dp), intent(in) :: y(:)
        integer, intent(in) :: ntrain
        type(garch_params_t), intent(in) :: params
        real(dp), intent(out) :: h(:)
        real(dp) :: backcast
        integer :: t

        backcast = sample_variance_n(y, ntrain)
        do t = 1, size(y)
            h(t) = max(params%omega + params%alpha*harch_block_local(y, t, 1, backcast) + &
                       params%gamma*harch_block_local(y, t, 5, backcast) + &
                       params%beta*harch_block_local(y, t, 22, backcast), 1.0e-12_dp)
        end do
    end subroutine harch_variance_path

    subroutine avgarch_variance_path(y, ntrain, params, h)
        real(dp), intent(in) :: y(:)
        integer, intent(in) :: ntrain
        type(garch_params_t), intent(in) :: params
        real(dp), intent(out) :: h(:)
        real(dp) :: sigma, z, x, news
        integer :: t

        sigma = max(params%omega / max(1.0_dp - avgarch_persist(params), 1.0e-8_dp), &
                    sqrt(sample_variance_n(y, ntrain)))
        do t = 1, size(y)
            sigma = max(sigma, 1.0e-6_dp)
            h(t) = sigma**2
            z = y(t) / sigma
            x = z - params%theta
            news = sqrt(1.0e-6_dp + x**2) - params%gamma*x
            sigma = params%omega + params%alpha*news*sigma + params%beta*sigma
        end do
    end subroutine avgarch_variance_path

    subroutine aparch_variance_path(y, ntrain, params, h)
        real(dp), intent(in) :: y(:)
        integer, intent(in) :: ntrain
        type(garch_params_t), intent(in) :: params
        real(dp), intent(out) :: h(:)
        real(dp) :: sdel, term
        integer :: t

        sdel = max(sample_variance_n(y, ntrain)**(0.5_dp*params%theta), 1.0e-12_dp)
        do t = 1, size(y)
            h(t) = max(sdel**(2.0_dp / params%theta), 1.0e-12_dp)
            term = max(abs(y(t)) - params%gamma*y(t), 1.0e-12_dp)**params%theta
            sdel = max(params%omega + params%alpha*term + params%beta*sdel, 1.0e-12_dp)
        end do
    end subroutine aparch_variance_path

    subroutine midas_variance_path(y, ntrain, params, h)
        real(dp), intent(in) :: y(:)
        integer, intent(in) :: ntrain
        type(garch_params_t), intent(in) :: params
        real(dp), intent(out) :: h(:)
        real(dp) :: weights(22), backcast
        integer :: t

        call midas_weights_local(params%theta, weights)
        backcast = sample_variance_n(y, ntrain)
        do t = 1, size(y)
            h(t) = max(params%omega + params%alpha*midas_block_local(y, t, weights, backcast), 1.0e-12_dp)
        end do
    end subroutine midas_variance_path

    subroutine midas_asym_variance_path(y, ntrain, params, h)
        real(dp), intent(in) :: y(:)
        integer, intent(in) :: ntrain
        type(garch_params_t), intent(in) :: params
        real(dp), intent(out) :: h(:)
        real(dp) :: weights(22), backcast, x, xneg
        integer :: t

        call midas_weights_local(params%theta, weights)
        backcast = sample_variance_n(y, ntrain)
        do t = 1, size(y)
            call midas_asym_blocks_local(y, t, weights, backcast, x, xneg)
            h(t) = max(params%omega + params%alpha*x + params%gamma*xneg, 1.0e-12_dp)
        end do
    end subroutine midas_asym_variance_path

    subroutine fgarch_twist_variance_path(y, params, h)
        real(dp), intent(in) :: y(:)
        type(garch_params_t), intent(in) :: params
        real(dp), intent(out) :: h(:)
        real(dp) :: moment, z, x, q
        integer :: t

        moment = fgarch_twist_moment_local(params%theta, params%twist)
        h(1) = max(params%omega / max(1.0_dp - params%alpha*moment - params%beta, 1.0e-8_dp), 1.0e-12_dp)
        do t = 1, size(y)
            h(t) = max(h(t), 1.0e-12_dp)
            z = y(t) / sqrt(h(t))
            x = z - params%theta
            q = abs(x) - params%twist*x
            if (t < size(y)) h(t + 1) = params%omega + params%alpha*h(t)*q**2 + params%beta*h(t)
        end do
    end subroutine fgarch_twist_variance_path

    real(dp) function sample_variance_n(y, n)
        real(dp), intent(in) :: y(:)
        integer, intent(in) :: n
        sample_variance_n = max(sum(y(1:n)**2) / real(n, dp), 1.0e-12_dp)
    end function sample_variance_n

    subroutine figarch_weights_local(phi, d, beta, lambda)
        real(dp), intent(in) :: phi, d, beta
        real(dp), intent(out) :: lambda(:)
        real(dp) :: delta_prev, delta_cur
        integer :: i

        delta_prev = d
        lambda(1) = d - phi + beta
        do i = 2, size(lambda)
            delta_cur = ((real(i, dp) - 1.0_dp - d) / real(i, dp))*delta_prev
            lambda(i) = beta*lambda(i - 1) + delta_cur - phi*delta_prev
            delta_prev = delta_cur
        end do
    end subroutine figarch_weights_local

    real(dp) function harch_block_local(y, t, lag, backcast)
        real(dp), intent(in) :: y(:), backcast
        integer, intent(in) :: t, lag
        integer :: j, idx

        harch_block_local = 0.0_dp
        do j = 1, lag
            idx = t - j
            if (idx >= 1) then
                harch_block_local = harch_block_local + y(idx)**2
            else
                harch_block_local = harch_block_local + backcast
            end if
        end do
        harch_block_local = harch_block_local / real(lag, dp)
    end function harch_block_local

    subroutine midas_weights_local(theta, weights)
        real(dp), intent(in) :: theta
        real(dp), intent(out) :: weights(:)
        real(dp) :: raw(size(weights)), ratio
        integer :: i

        raw(1) = theta
        do i = 2, size(weights)
            ratio = (real(i - 1, dp) + theta) / real(i, dp)
            raw(i) = raw(i - 1)*ratio
        end do
        weights = raw / max(sum(raw), 1.0e-12_dp)
    end subroutine midas_weights_local

    real(dp) function midas_block_local(y, t, weights, backcast)
        real(dp), intent(in) :: y(:), weights(:), backcast
        integer, intent(in) :: t
        integer :: i, idx

        midas_block_local = 0.0_dp
        do i = 1, size(weights)
            idx = t - i
            if (idx >= 1) then
                midas_block_local = midas_block_local + weights(i)*y(idx)**2
            else
                midas_block_local = midas_block_local + weights(i)*backcast
            end if
        end do
    end function midas_block_local

    subroutine midas_asym_blocks_local(y, t, weights, backcast, x, xneg)
        real(dp), intent(in) :: y(:), weights(:), backcast
        integer, intent(in) :: t
        real(dp), intent(out) :: x, xneg
        real(dp) :: lag_sq, ind
        integer :: i, idx

        x = 0.0_dp
        xneg = 0.0_dp
        do i = 1, size(weights)
            idx = t - i
            if (idx >= 1) then
                lag_sq = y(idx)**2
                ind = merge(1.0_dp, 0.0_dp, y(idx) < 0.0_dp)
            else
                lag_sq = backcast
                ind = 0.5_dp
            end if
            x = x + weights(i)*lag_sq
            xneg = xneg + weights(i)*ind*lag_sq
        end do
    end subroutine midas_asym_blocks_local

    real(dp) function fgarch_twist_moment_local(theta, c)
        real(dp), intent(in) :: theta, c
        real(dp) :: phi, ph, a, b

        phi = 0.5_dp * (1.0_dp + erf(theta / sqrt(2.0_dp)))
        ph = exp(-0.5_dp*theta**2) / sqrt(2.0_dp*acos(-1.0_dp))
        a = (1.0_dp + theta**2)*(1.0_dp - phi) - theta*ph
        b = (1.0_dp + theta**2)*phi + theta*ph
        fgarch_twist_moment_local = (1.0_dp - c)**2*a + (1.0_dp + c)**2*b
    end function fgarch_twist_moment_local

    subroutine set_true_params(model_name, params)
        character(len=*), intent(in) :: model_name
        type(garch_params_t), intent(out) :: params

        params = garch_params_t()

        select case (trim(adjustl(model_name)))
        case ("SYMM_GARCH")
            params%omega = 1.0e-6_dp
            params%alpha = 0.0700_dp
            params%beta = 0.9000_dp
        case ("SYMM_GARCH21")
            params%omega = 1.0e-6_dp
            allocate(params%alpha_lags(2), params%beta_lags(1))
            params%alpha_lags = [0.0500_dp, 0.0200_dp]
            params%beta_lags = [0.9000_dp]
            params%alpha = sum(params%alpha_lags)
            params%beta = sum(params%beta_lags)
        case ("NAGARCH")
            params%omega = 1.2e-6_dp
            params%alpha = 0.0700_dp
            params%beta = 0.8800_dp
            params%theta = 0.4000_dp
        case ("NAGARCH21")
            params%omega = 1.2e-6_dp
            allocate(params%alpha_lags(2), params%beta_lags(1))
            params%alpha_lags = [0.0400_dp, 0.0200_dp]
            params%beta_lags = [0.8600_dp]
            params%alpha = sum(params%alpha_lags)
            params%beta = sum(params%beta_lags)
            params%theta = 0.3500_dp
        case ("GJR_GARCH")
            params%omega = 1.0e-6_dp
            params%alpha = 0.0500_dp
            params%gamma = 0.0800_dp
            params%beta = 0.8800_dp
        case ("GJR_SIGNED")
            params%omega = 1.0e-6_dp
            params%alpha = 0.0800_dp
            params%gamma = -0.0400_dp
            params%beta = 0.8900_dp
        case ("EGARCH")
            params%omega = -0.1500_dp
            params%alpha = 0.1200_dp
            params%gamma = -0.0500_dp
            params%beta = 0.9500_dp
        case ("QGARCH")
            params%omega = 1.0e-6_dp
            params%alpha = 0.0800_dp
            params%beta = 0.8600_dp
            params%theta = -0.0400_dp
        case ("CSGARCH")
            params%omega = 1.0e-6_dp
            params%alpha = 0.0600_dp
            params%beta = 0.7500_dp
            params%gamma = 0.0300_dp
            params%extra1 = 0.9700_dp
            params%extra2 = 0.0300_dp
        case ("TGARCH")
            params%omega = 1.0e-6_dp
            params%alpha = 0.0600_dp
            params%beta = 0.9000_dp
            params%gamma = -0.2500_dp
        case ("FIGARCH")
            params%omega = 1.0e-6_dp
            params%alpha = 0.1200_dp
            params%theta = 0.3500_dp
            params%beta = 0.3000_dp
        case ("FI_NAGARCH")
            params%omega = 1.0e-6_dp
            params%alpha = 0.1000_dp
            params%theta = 0.3000_dp
            params%beta = 0.3000_dp
            params%twist = 0.3000_dp
        case ("HARCH")
            params%omega = 1.0e-6_dp
            params%alpha = 0.0800_dp
            params%gamma = 0.3000_dp
            params%beta = 0.5500_dp
        case ("AVGARCH")
            params%omega = 1.0e-6_dp
            params%alpha = 0.0600_dp
            params%beta = 0.9000_dp
            params%gamma = -0.2000_dp
            params%theta = 0.0500_dp
        case ("APARCH")
            params%omega = 1.0e-6_dp
            params%alpha = 0.0700_dp
            params%gamma = 0.4000_dp
            params%beta = 0.8800_dp
            params%theta = 1.2000_dp
        case ("MIDAS_HYPER")
            params%omega = 1.0e-6_dp
            params%alpha = 0.9000_dp
            params%theta = 0.5000_dp
        case ("MIDAS_ASYM")
            params%omega = 1.0e-6_dp
            params%alpha = 0.7500_dp
            params%gamma = 0.2500_dp
            params%theta = 0.5000_dp
        case ("FGARCH_TWIST")
            params%omega = 1.0e-6_dp
            params%alpha = 0.0600_dp
            params%beta = 0.8800_dp
            params%theta = 0.1000_dp
            params%twist = -0.2000_dp
        case default
            error stop "xsim_garch_fit: no true parameters for model"
        end select
    end subroutine set_true_params

    subroutine fit_one_model(model_name, y, max_iter, gtol, f_best, params, niter, converged)
        character(len=*), intent(in) :: model_name
        real(dp), intent(in) :: y(:), gtol
        integer, intent(in) :: max_iter
        real(dp), intent(out) :: f_best
        type(garch_params_t), intent(out) :: params
        integer, intent(out) :: niter
        logical, intent(out) :: converged
        real(dp) :: vol_ann, skew, ekurt

        select case (trim(adjustl(model_name)))
        case ("SYMM_GARCH")
            call fit_symm_garch(y, max_iter, gtol, f_best, params, niter, converged)
        case ("SYMM_GARCH21")
            call fit_symm_garch_pq(y, 2, 1, max_iter, gtol, f_best, params, niter, converged)
        case ("NAGARCH")
            call fit_nagarch(y, max_iter, gtol, f_best, params, niter, converged)
        case ("NAGARCH21")
            call fit_nagarch_pq(y, 2, 1, max_iter, gtol, f_best, params, niter, converged)
        case ("GJR_GARCH")
            call fit_gjr(y, max_iter, gtol, f_best, params, niter, converged)
        case ("GJR_SIGNED")
            call fit_gjr_signed(y, max_iter, gtol, f_best, params, niter, converged)
        case ("EGARCH")
            call fit_egarch(y, max_iter, gtol, f_best, params, niter, converged)
        case ("QGARCH")
            call fit_qgarch(y, max_iter, gtol, f_best, params, niter, converged)
        case ("CSGARCH")
            call fit_csgarch(y, max_iter, gtol, f_best, params, niter, converged)
        case ("TGARCH")
            call fit_tgarch(y, max_iter, gtol, f_best, params, niter, converged)
        case ("FIGARCH")
            call fit_figarch(y, max_iter, gtol, f_best, params, niter, converged)
        case ("FI_NAGARCH")
            call fit_fi_nagarch(y, max_iter, gtol, f_best, params, niter, converged)
        case ("HARCH")
            call fit_harch(y, max_iter, gtol, f_best, params, niter, converged)
        case ("AVGARCH")
            call fit_avgarch(y, max_iter, gtol, f_best, params, niter, converged)
        case ("APARCH")
            call fit_aparch(y, max_iter, gtol, f_best, params, niter, converged)
        case ("MIDAS_HYPER")
            call fit_midas_hyperbolic(y, max_iter, gtol, f_best, params, niter, converged)
        case ("MIDAS_ASYM")
            call fit_midas_hyperbolic_asym(y, max_iter, gtol, f_best, params, niter, converged)
        case ("FGARCH_TWIST")
            call fit_fgarch_twist(y, sqrt(max(sum(y**2) / real(size(y), dp), 1.0e-12_dp)), &
                                  max_iter, gtol, f_best, params, vol_ann, skew, ekurt, &
                                  niter, converged)
        case default
            error stop "xsim_garch_fit: no fitter for model"
        end select
    end subroutine fit_one_model

    real(dp) function model_persist(model_name, params)
        character(len=*), intent(in) :: model_name
        type(garch_params_t), intent(in) :: params

        select case (trim(adjustl(model_name)))
        case ("SYMM_GARCH")
            model_persist = symm_garch_persist(params)
        case ("SYMM_GARCH21")
            model_persist = symm_garch_pq_persist(params)
        case ("NAGARCH")
            model_persist = nagarch_persist(params)
        case ("NAGARCH21")
            model_persist = nagarch_pq_persist(params)
        case ("GJR_GARCH", "GJR_SIGNED")
            model_persist = gjr_persist(params)
        case ("EGARCH")
            model_persist = egarch_persist(params)
        case ("QGARCH")
            model_persist = qgarch_persist(params)
        case ("CSGARCH")
            model_persist = csgarch_persist(params)
        case ("TGARCH")
            model_persist = tgarch_persist(params)
        case ("FIGARCH")
            model_persist = figarch_persist(params)
        case ("FI_NAGARCH")
            model_persist = fi_nagarch_persist(params)
        case ("HARCH")
            model_persist = harch_persist(params)
        case ("AVGARCH")
            model_persist = avgarch_persist(params)
        case ("APARCH")
            model_persist = aparch_persist(params)
        case ("MIDAS_HYPER")
            model_persist = midas_hyperbolic_persist(params)
        case ("MIDAS_ASYM")
            model_persist = midas_hyperbolic_asym_persist(params)
        case ("FGARCH_TWIST")
            model_persist = fgarch_twist_persist(params)
        case default
            error stop "xsim_garch_fit: no persistence function for model"
        end select
    end function model_persist

end program xsim_garch_fit

program xsim_garch_fit_dist
    use, intrinsic :: iso_fortran_env, only: int64
    use kind_mod, only: dp
    use stats_mod, only: variance
    use distributions_mod, only: pdf_normal, pdf_t, pdf_ged, pdf_logistic, &
                                 pdf_laplace, pdf_sech, pdf_nig, pdf_fs_skewt
    use random_mod, only: random_normal, random_t_std, random_ged_std, &
                          random_laplace_std, random_logistic_std, random_sech, &
                          random_nig_sym, random_fs_skewt
    use garch_types_mod, only: garch_params_t
    use garch_fit_dist_mod, only: garch_dist_fit_result_t, fit_garch_dist_model
    implicit none

    integer, parameter :: dist_normal = 1
    integer, parameter :: dist_t = 2
    integer, parameter :: dist_sech = 3
    integer, parameter :: dist_ged = 4
    integer, parameter :: dist_laplace = 5
    integer, parameter :: dist_logistic = 6
    integer, parameter :: dist_nig = 7
    integer, parameter :: dist_fs_skewt = 8
    integer, parameter :: n_model = 6
    integer, parameter :: n_dist = 8
    integer, parameter :: n_keep = 500
    integer, parameter :: n_burn = 100
    integer, parameter :: max_iter = 500
    real(dp), parameter :: gtol = 1.0e-7_dp
    real(dp), parameter :: true_omega = 1.0e-6_dp
    real(dp), parameter :: true_alpha = 0.07_dp
    real(dp), parameter :: true_beta = 0.90_dp
    real(dp), parameter :: min_pdf = 1.0e-300_dp

    type :: fit_result_t
        real(dp) :: omega = 0.0_dp
        real(dp) :: alpha = 0.0_dp
        real(dp) :: gamma = 0.0_dp
        real(dp) :: beta = 0.0_dp
        real(dp) :: theta = 0.0_dp
        real(dp) :: shape = 0.0_dp
        real(dp) :: shape2 = 0.0_dp
        real(dp) :: nll = huge(1.0_dp)
        integer :: niter = 0
        logical :: converged = .false.
    end type fit_result_t

    character(len=10), parameter :: dist_names(n_dist) = [ &
        "NORMAL    ", "T         ", "SECH      ", "GED       ", "LAPLACE   ", "LOGISTIC  ", "NIG       ", "FS_SKEWT  " ]
    character(len=16), parameter :: model_names(n_model) = [ character(len=16) :: &
        "SYMM_GARCH", "NAGARCH", "GJR_GARCH", "GJR_SIGNED", "EGARCH", "QGARCH" ]
    real(dp), parameter :: true_shape(n_dist) = [ &
        0.0_dp, 6.0_dp, 0.0_dp, 1.5_dp, 0.0_dp, 0.0_dp, 3.0_dp, 6.0_dp ]
    real(dp), parameter :: true_shape2(n_dist) = [ &
        0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 1.2_dp ]

    real(dp), allocatable :: y_all(:), y(:), y_train_test(:)
    type(fit_result_t) :: fit_true, fit_all(n_dist)
    type(fit_result_t) :: true_fit
    real(dp) :: true_loss, fit_loss(n_dist)
    real(dp) :: win_vs_normal(n_model,n_dist), win_vs_best_other(n_model,n_dist)
    real(dp) :: sim_time(n_model,n_dist), fit_time(n_model,n_dist)
    real(dp) :: sum_sim_time, sum_fit_time
    integer :: imodel, idist, jdist, irep, n_rep, ntest, n_total, ios, seed_val
    integer :: true_rank, best_other_idx
    integer(int64) :: clock_count0, clock_count1, clock_rate
    character(len=64) :: arg1

    n_rep = 1
    ntest = 0
    call system_clock(count_rate=clock_rate)

    if (command_argument_count() >= 1) then
        call get_command_argument(1, arg1)
        read(arg1, *, iostat=ios) n_rep
        if (ios /= 0 .or. n_rep < 1) then
            write(*,'(A)') "Usage: xsim_garch_fit_dist.exe [n_rep>=1] [ntest>=0]"
            error stop "xsim_garch_fit_dist: invalid n_rep argument"
        end if
    end if
    if (command_argument_count() >= 2) then
        call get_command_argument(2, arg1)
        read(arg1, *, iostat=ios) ntest
        if (ios /= 0 .or. ntest < 0) then
            write(*,'(A)') "Usage: xsim_garch_fit_dist.exe [n_rep>=1] [ntest>=0]"
            error stop "xsim_garch_fit_dist: invalid ntest argument"
        end if
    end if

    n_total = n_keep + n_burn + ntest
    allocate(y_all(n_total), y(n_keep), y_train_test(max(1, n_keep + ntest)))

    win_vs_normal = 0.0_dp
    win_vs_best_other = 0.0_dp
    sim_time = 0.0_dp
    fit_time = 0.0_dp

    write(*,'(A)') "GARCH(1,1) distribution simulate -> fit sanity check"
    write(*,'(A,I0)') "Observations simulated per dataset (incl burn-in): ", n_total
    write(*,'(A,I0,A,I0,A,I0,A,I0,A,I0)') "Kept observations: ", n_keep, &
        "  burn-in: ", n_burn, "  datasets/model-dist: ", n_rep, &
        "  models: ", n_model, "  distributions: ", n_dist
    if (ntest > 0) then
        write(*,'(A,I0,A)') "Out-of-sample observations: ", ntest, " (enabled)"
    else
        write(*,'(A,I0,A)') "Out-of-sample observations: ", ntest, " (disabled)"
    end if
    write(*,'(A)') "Rows: TRUE parameters and FIT estimates from simulated returns"
    write(*,'(A)') repeat("-", 162)
    write(*,'(A)') "Model            Dist       Row   Set      omega    alpha    gamma     beta    theta      nu      xi  dist_alpha  persist  niter conv    nll/obs"
    write(*,'(A)') repeat("-", 162)

    do imodel = 1, n_model
        do idist = 1, n_dist
            sum_sim_time = 0.0_dp
            sum_fit_time = 0.0_dp
            call set_true_model(trim(model_names(imodel)), true_fit)
            true_fit%shape = true_shape(idist)
            true_fit%shape2 = true_shape2(idist)

            write(*,'(A16,1X,A10,1X,A5,1X,A4,1X,ES9.2,1X,F8.4,1X,F8.4,1X,F8.4,1X,F8.4,1X,A8,1X,A8,1X,A10,1X,F8.4,1X,A5,1X,A4,1X,A10)') &
                trim(model_names(imodel)), trim(dist_names(idist)), "TRUE", "-", true_fit%omega, &
                true_fit%alpha, true_fit%gamma, true_fit%beta, true_fit%theta, &
                dist_nu_value(idist, true_fit%shape), dist_xi_value(idist, true_fit%shape2), &
                dist_alpha_value(idist, true_fit%shape), &
                model_persist(trim(model_names(imodel)), true_fit), "-", "-", "-"

            do irep = 1, n_rep
                seed_val = 200000 + 10000*imodel + 1000*idist + irep

                call system_clock(clock_count0)
                call simulate_dist_model(trim(model_names(imodel)), idist, true_fit, n_total, seed_val, y_all)
                call system_clock(clock_count1)
                sum_sim_time = sum_sim_time + elapsed_seconds(clock_count0, clock_count1, clock_rate)

                y = y_all(n_burn + 1:n_burn + n_keep)
                if (ntest > 0) y_train_test = y_all(n_burn + 1:n_total)

                call system_clock(clock_count0)
                call fit_dist_model(trim(model_names(imodel)), idist, y, n_keep, max_iter, gtol, fit_true)
                call system_clock(clock_count1)
                sum_fit_time = sum_fit_time + elapsed_seconds(clock_count0, clock_count1, clock_rate)

                write(*,'(A16,1X,A10,1X,A5,1X,I4,1X,ES9.2,1X,F8.4,1X,F8.4,1X,F8.4,1X,F8.4,1X,A8,1X,A8,1X,A10,1X,F8.4,1X,I5,1X,L4,1X,F10.6)') &
                    trim(model_names(imodel)), trim(dist_names(idist)), "FIT", irep, fit_true%omega, &
                    fit_true%alpha, fit_true%gamma, fit_true%beta, fit_true%theta, &
                    dist_nu_value(idist, fit_true%shape), dist_xi_value(idist, fit_true%shape2), &
                    dist_alpha_value(idist, fit_true%shape), &
                    model_persist(trim(model_names(imodel)), fit_true), &
                    fit_true%niter, fit_true%converged, fit_true%nll

                if (ntest > 0) then
                    do jdist = 1, n_dist
                        if (jdist == idist) then
                            fit_all(jdist) = fit_true
                        else
                            call fit_dist_model(trim(model_names(imodel)), jdist, y, n_keep, max_iter, gtol, fit_all(jdist))
                        end if
                        fit_loss(jdist) = oos_nll_dist(trim(model_names(imodel)), jdist, y_train_test, n_keep, &
                                                       ntest, fit_all(jdist))
                    end do
                    true_loss = fit_loss(idist)
                    if (true_loss < fit_loss(dist_normal)) then
                        win_vs_normal(imodel,idist) = win_vs_normal(imodel,idist) + 1.0_dp
                    end if
                    best_other_idx = best_other_dist(idist, fit_loss)
                    if (true_loss < fit_loss(best_other_idx)) then
                        win_vs_best_other(imodel,idist) = win_vs_best_other(imodel,idist) + 1.0_dp
                    end if
                    true_rank = 1 + count(fit_loss < true_loss)
                    call print_oos_row(trim(model_names(imodel)), idist, fit_loss, true_rank)
                end if
            end do

            sim_time(imodel,idist) = sum_sim_time / real(n_rep, dp)
            fit_time(imodel,idist) = sum_fit_time / real(n_rep, dp)
            write(*,'(A)') repeat("-", 162)
        end do
    end do

    if (ntest > 0) then
        call print_oos_win_table(win_vs_normal / real(n_rep, dp), &
                                 win_vs_best_other / real(n_rep, dp), sim_time, fit_time)
    end if

contains

    real(dp) function elapsed_seconds(count0, count1, count_rate)
        integer(int64), intent(in) :: count0, count1, count_rate

        if (count_rate <= 0_int64 .or. count1 < count0) then
            elapsed_seconds = 0.0_dp
        else
            elapsed_seconds = real(count1 - count0, dp) / real(count_rate, dp)
        end if
    end function elapsed_seconds

    subroutine seed_rng(seed_val)
        integer, intent(in) :: seed_val
        integer :: sz
        integer, allocatable :: seed_arr(:)

        call random_seed(size=sz)
        allocate(seed_arr(sz))
        seed_arr = seed_val
        call random_seed(put=seed_arr)
        deallocate(seed_arr)
    end subroutine seed_rng

    character(len=8) function dist_nu_value(dist_id, shape)
        integer, intent(in) :: dist_id
        real(dp), intent(in) :: shape

        select case (dist_id)
        case (dist_fs_skewt)
            write(dist_nu_value, '(F8.3)') shape
        case (dist_t, dist_ged)
            write(dist_nu_value, '(F8.3)') shape
        case default
            dist_nu_value = "-"
        end select
    end function dist_nu_value

    character(len=8) function dist_xi_value(dist_id, shape2)
        integer, intent(in) :: dist_id
        real(dp), intent(in) :: shape2

        select case (dist_id)
        case (dist_fs_skewt)
            write(dist_xi_value, '(F8.3)') shape2
        case default
            dist_xi_value = "-"
        end select
    end function dist_xi_value

    character(len=10) function dist_alpha_value(dist_id, shape)
        integer, intent(in) :: dist_id
        real(dp), intent(in) :: shape

        select case (dist_id)
        case (dist_nig)
            write(dist_alpha_value, '(F10.3)') shape
        case default
            dist_alpha_value = "-"
        end select
    end function dist_alpha_value

    real(dp) function draw_innovation(dist_id, shape, shape2)
        integer, intent(in) :: dist_id
        real(dp), intent(in) :: shape, shape2

        select case (dist_id)
        case (dist_normal)
            draw_innovation = random_normal()
        case (dist_t)
            draw_innovation = random_t_std(shape)
        case (dist_sech)
            draw_innovation = random_sech()
        case (dist_ged)
            draw_innovation = random_ged_std(shape)
        case (dist_laplace)
            draw_innovation = random_laplace_std()
        case (dist_logistic)
            draw_innovation = random_logistic_std()
        case (dist_nig)
            draw_innovation = random_nig_sym(shape)
        case (dist_fs_skewt)
            draw_innovation = random_fs_skewt(shape, shape2)
        case default
            error stop "draw_innovation: unknown distribution"
        end select
    end function draw_innovation

    subroutine set_true_model(model_name, params)
        character(len=*), intent(in) :: model_name
        type(fit_result_t), intent(out) :: params

        params = fit_result_t()
        select case (model_name)
        case ("SYMM_GARCH", "GARCH")
            params%omega = true_omega
            params%alpha = true_alpha
            params%beta = true_beta
        case ("NAGARCH")
            params%omega = 1.2e-6_dp
            params%alpha = 0.07_dp
            params%beta = 0.88_dp
            params%theta = 0.40_dp
        case ("GJR_GARCH", "GJR")
            params%omega = 1.0e-6_dp
            params%alpha = 0.05_dp
            params%gamma = 0.08_dp
            params%beta = 0.88_dp
        case ("GJR_SIGNED")
            params%omega = 1.0e-6_dp
            params%alpha = 0.08_dp
            params%gamma = -0.04_dp
            params%beta = 0.89_dp
        case ("EGARCH")
            params%omega = -0.15_dp
            params%alpha = 0.12_dp
            params%gamma = -0.05_dp
            params%beta = 0.95_dp
        case ("QGARCH")
            params%omega = 1.0e-6_dp
            params%alpha = 0.08_dp
            params%beta = 0.86_dp
            params%theta = -0.04_dp
        case default
            error stop "set_true_model: unsupported model"
        end select
    end subroutine set_true_model

    real(dp) function model_persist(model_name, params)
        character(len=*), intent(in) :: model_name
        type(fit_result_t), intent(in) :: params

        select case (model_name)
        case ("SYMM_GARCH", "GARCH")
            model_persist = params%alpha + params%beta
        case ("NAGARCH")
            model_persist = params%beta + params%alpha*(1.0_dp + params%theta**2)
        case ("GJR_GARCH", "GJR")
            model_persist = params%alpha + 0.5_dp*params%gamma + params%beta
        case ("GJR_SIGNED")
            model_persist = params%alpha + params%beta
        case ("EGARCH")
            model_persist = params%beta
        case ("QGARCH")
            model_persist = params%alpha + params%beta
        case default
            model_persist = params%alpha + params%beta
        end select
    end function model_persist

    subroutine simulate_dist_model(model_name, dist_id, params, n, seed_val, y)
        character(len=*), intent(in) :: model_name
        integer, intent(in) :: dist_id, n, seed_val
        type(fit_result_t), intent(in) :: params
        real(dp), intent(out) :: y(n)
        real(dp) :: h, lh, eps, sqrth, ind, c_eg
        integer :: i

        select case (model_name)
        case ("SYMM_GARCH", "GARCH")
            call seed_rng(seed_val)
            h = max(params%omega / max(1.0_dp - params%alpha - params%beta, 1.0e-8_dp), 1.0e-12_dp)
            do i = 1, n
                eps = draw_innovation(dist_id, params%shape, params%shape2)
                y(i) = sqrt(h) * eps
                h = max(params%omega + params%alpha*y(i)**2 + params%beta*h, 1.0e-12_dp)
            end do
        case ("NAGARCH")
            call seed_rng(seed_val)
            h = max(params%omega / max(1.0_dp - model_persist(model_name, params), 1.0e-8_dp), 1.0e-12_dp)
            do i = 1, n
                eps = draw_innovation(dist_id, params%shape, params%shape2)
                sqrth = sqrt(h)
                y(i) = sqrth * eps
                h = max(params%omega + params%alpha*(y(i) - params%theta*sqrth)**2 + params%beta*h, 1.0e-12_dp)
            end do
        case ("GJR_GARCH", "GJR")
            call seed_rng(seed_val)
            h = max(params%omega / max(1.0_dp - model_persist(model_name, params), 1.0e-8_dp), 1.0e-12_dp)
            do i = 1, n
                eps = draw_innovation(dist_id, params%shape, params%shape2)
                y(i) = sqrt(h) * eps
                ind = merge(1.0_dp, 0.0_dp, y(i) < 0.0_dp)
                h = max(params%omega + (params%alpha + params%gamma*ind)*y(i)**2 + params%beta*h, 1.0e-12_dp)
            end do
        case ("GJR_SIGNED")
            call seed_rng(seed_val)
            h = max(params%omega / max(1.0_dp - model_persist(model_name, params), 1.0e-8_dp), 1.0e-12_dp)
            do i = 1, n
                eps = draw_innovation(dist_id, params%shape, params%shape2)
                y(i) = sqrt(h) * eps
                ind = sign(1.0_dp, y(i))
                h = max(params%omega + max(params%alpha + params%gamma*ind, 0.0_dp)*y(i)**2 + &
                        params%beta*h, 1.0e-12_dp)
            end do
        case ("EGARCH")
            call seed_rng(seed_val)
            c_eg = sqrt(2.0_dp / acos(-1.0_dp))
            lh = params%omega / max(1.0_dp - params%beta, 1.0e-8_dp)
            do i = 1, n
                h = max(exp(lh), 1.0e-12_dp)
                eps = draw_innovation(dist_id, params%shape, params%shape2)
                y(i) = sqrt(h) * eps
                lh = params%omega + params%beta*lh + params%alpha*(abs(eps) - c_eg) + params%gamma*eps
            end do
        case ("QGARCH")
            call seed_rng(seed_val)
            h = max((params%omega + params%alpha*params%theta**2) / &
                    max(1.0_dp - params%alpha - params%beta, 1.0e-8_dp), 1.0e-12_dp)
            do i = 1, n
                eps = draw_innovation(dist_id, params%shape, params%shape2)
                y(i) = sqrt(h) * eps
                h = max(params%omega + params%alpha*(y(i) - params%theta)**2 + params%beta*h, 1.0e-12_dp)
            end do
        case default
            error stop "simulate_dist_model: unsupported model"
        end select
    end subroutine simulate_dist_model

    subroutine fit_dist_model(model_name, dist_id, y, n, max_iter, gtol, fit)
        character(len=*), intent(in) :: model_name
        integer, intent(in) :: dist_id, n, max_iter
        real(dp), intent(in) :: y(n), gtol
        type(fit_result_t), intent(out) :: fit
        type(garch_dist_fit_result_t) :: result

        call fit_garch_dist_model(model_name, trim(dist_names(dist_id)), y, max_iter, gtol, result)
        fit%omega = result%params%omega
        fit%alpha = result%params%alpha
        fit%gamma = result%params%gamma
        fit%beta = result%params%beta
        fit%theta = result%params%theta
        fit%shape = result%shape
        fit%shape2 = result%shape2
        fit%nll = result%nll
        fit%niter = result%niter
        fit%converged = result%converged
    end subroutine fit_dist_model

    subroutine variance_path(model_name, y, fit, h)
        character(len=*), intent(in) :: model_name
        real(dp), intent(in) :: y(:)
        type(fit_result_t), intent(in) :: fit
        real(dp), intent(out) :: h(:)
        real(dp) :: lh, z, c_eg, sqrth, ind
        integer :: t

        select case (model_name)
        case ("SYMM_GARCH", "GARCH")
            h(1) = max(fit%omega / max(1.0_dp - fit%alpha - fit%beta, 1.0e-8_dp), variance(y))
            h(1) = max(h(1), 1.0e-12_dp)
            do t = 2, size(y)
                h(t) = max(fit%omega + fit%alpha*y(t - 1)**2 + fit%beta*h(t - 1), 1.0e-12_dp)
            end do
        case ("NAGARCH")
            h(1) = max(fit%omega / max(1.0_dp - model_persist(model_name, fit), 1.0e-8_dp), variance(y))
            h(1) = max(h(1), 1.0e-12_dp)
            do t = 2, size(y)
                sqrth = sqrt(max(h(t - 1), 1.0e-12_dp))
                h(t) = max(fit%omega + fit%alpha*(y(t - 1) - fit%theta*sqrth)**2 + &
                           fit%beta*h(t - 1), 1.0e-12_dp)
            end do
        case ("GJR_GARCH", "GJR")
            h(1) = max(fit%omega / max(1.0_dp - model_persist(model_name, fit), 1.0e-8_dp), variance(y))
            h(1) = max(h(1), 1.0e-12_dp)
            do t = 2, size(y)
                ind = merge(1.0_dp, 0.0_dp, y(t - 1) < 0.0_dp)
                h(t) = max(fit%omega + (fit%alpha + fit%gamma*ind)*y(t - 1)**2 + &
                           fit%beta*h(t - 1), 1.0e-12_dp)
            end do
        case ("GJR_SIGNED")
            h(1) = max(fit%omega / max(1.0_dp - model_persist(model_name, fit), 1.0e-8_dp), variance(y))
            h(1) = max(h(1), 1.0e-12_dp)
            do t = 2, size(y)
                ind = sign(1.0_dp, y(t - 1))
                h(t) = max(fit%omega + max(fit%alpha + fit%gamma*ind, 0.0_dp)*y(t - 1)**2 + &
                           fit%beta*h(t - 1), 1.0e-12_dp)
            end do
        case ("EGARCH")
            c_eg = sqrt(2.0_dp / acos(-1.0_dp))
            lh = fit%omega / max(1.0_dp - fit%beta, 1.0e-8_dp)
            do t = 1, size(y)
                h(t) = max(exp(lh), 1.0e-12_dp)
                z = y(t) / sqrt(h(t))
                lh = fit%omega + fit%beta*lh + fit%alpha*(abs(z) - c_eg) + fit%gamma*z
            end do
        case ("QGARCH")
            h(1) = max((fit%omega + fit%alpha*fit%theta**2) / &
                       max(1.0_dp - fit%alpha - fit%beta, 1.0e-8_dp), variance(y))
            h(1) = max(h(1), 1.0e-12_dp)
            do t = 2, size(y)
                h(t) = max(fit%omega + fit%alpha*(y(t - 1) - fit%theta)**2 + &
                           fit%beta*h(t - 1), 1.0e-12_dp)
            end do
        case default
            error stop "variance_path: unsupported model"
        end select
    end subroutine variance_path

    real(dp) function oos_nll_dist(model_name, dist_id, y_full, ntrain, ntest, fit)
        character(len=*), intent(in) :: model_name
        integer, intent(in) :: dist_id, ntrain, ntest
        real(dp), intent(in) :: y_full(:)
        type(fit_result_t), intent(in) :: fit
        real(dp), allocatable :: h(:)
        real(dp) :: z, dens, loss
        integer :: t

        allocate(h(size(y_full)))
        call variance_path(model_name, y_full, fit, h)
        loss = 0.0_dp
        do t = ntrain + 1, ntrain + ntest
            z = y_full(t) / sqrt(max(h(t), 1.0e-12_dp))
            dens = max(innovation_pdf(dist_id, z, fit%shape, fit%shape2), min_pdf)
            loss = loss - log(dens) + 0.5_dp * log(max(h(t), 1.0e-12_dp))
        end do
        oos_nll_dist = loss / real(ntest, dp)
        deallocate(h)
    end function oos_nll_dist

    real(dp) function innovation_pdf(dist_id, z, shape, shape2)
        integer, intent(in) :: dist_id
        real(dp), intent(in) :: z, shape, shape2

        select case (dist_id)
        case (dist_normal)
            innovation_pdf = pdf_normal(z)
        case (dist_t)
            innovation_pdf = pdf_t(z, shape)
        case (dist_sech)
            innovation_pdf = pdf_sech(z)
        case (dist_ged)
            innovation_pdf = pdf_ged(z, shape)
        case (dist_laplace)
            innovation_pdf = pdf_laplace(z)
        case (dist_logistic)
            innovation_pdf = pdf_logistic(z)
        case (dist_nig)
            innovation_pdf = pdf_nig(z, shape)
        case (dist_fs_skewt)
            innovation_pdf = pdf_fs_skewt(z, shape, shape2)
        case default
            innovation_pdf = min_pdf
        end select
    end function innovation_pdf

    integer function best_other_dist(true_dist, loss)
        integer, intent(in) :: true_dist
        real(dp), intent(in) :: loss(n_dist)
        real(dp) :: best_loss
        integer :: j

        best_other_dist = 0
        best_loss = huge(1.0_dp)
        do j = 1, n_dist
            if (j == true_dist) cycle
            if (loss(j) < best_loss) then
                best_loss = loss(j)
                best_other_dist = j
            end if
        end do
    end function best_other_dist

    subroutine print_oos_row(model_name, true_dist, loss, true_rank)
        character(len=*), intent(in) :: model_name
        integer, intent(in) :: true_dist, true_rank
        real(dp), intent(in) :: loss(n_dist)
        integer :: j

        write(*,'(A,1X,A,1X,A,1X,A,I0)') "OOS_NLL", trim(model_name), trim(dist_names(true_dist)), &
            "true_fit_rank=", true_rank
        write(*,'(10X,A10,1X,A11)') "Fit_dist", "mean_nll"
        do j = 1, n_dist
            write(*,'(10X,A10,1X,F11.6)') dist_names(j), loss(j)
        end do
    end subroutine print_oos_row

    subroutine print_oos_win_table(win_normal, win_best_other, sim_time, fit_time)
        real(dp), intent(in) :: win_normal(n_model,n_dist), win_best_other(n_model,n_dist)
        real(dp), intent(in) :: sim_time(n_model,n_dist), fit_time(n_model,n_dist)
        integer :: im, j

        write(*,'(/,A)') "OOS true_fit_win_rate by simulated model and distribution"
        write(*,'(A)') repeat("-", 104)
        write(*,'(A16,1X,A10,1X,A12,1X,A18,1X,A12,1X,A12)') &
            "Model", "Dist", "vs_NORMAL", "vs_best_other", "sim_sec", "fit_sec"
        write(*,'(A)') repeat("-", 104)
        do im = 1, n_model
            do j = 1, n_dist
                write(*,'(A16,1X,A10,1X,F12.3,1X,F18.3,1X,F12.6,1X,F12.6)') &
                    trim(model_names(im)), trim(dist_names(j)), win_normal(im,j), &
                    win_best_other(im,j), sim_time(im,j), fit_time(im,j)
            end do
        end do
        write(*,'(A)') repeat("-", 104)
    end subroutine print_oos_win_table

end program xsim_garch_fit_dist

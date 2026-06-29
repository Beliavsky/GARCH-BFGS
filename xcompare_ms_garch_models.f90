! Full-sample AIC/BIC comparison of GARCH, NAGARCH, MSGARCH, and MSNAGARCH
! across one or more noise distributions (default: normal and t).
!
! For each distribution, all four models are fitted with that distribution and
! AIC/BIC are computed from the fitted log-likelihood, giving an apples-to-apples
! comparison within each table.
!
! Usage: xcompare_ms_garch_models [prices_csv] [-numerical] [-dist dist1[,dist2...]]
!   dist names: NORMAL, T  (case-insensitive; default: NORMAL,T)
!   Example: xcompare_ms_garch_models vix_spy.csv -dist T

program xcompare_ms_garch_models
    use kind_mod,             only: dp
    use date_mod,             only: print_program_header
    use csv_mod,              only: read_price_csv, print_price_sample_info
    use stats_mod,            only: mean
    use garch_types_mod,      only: garch_params_t
    use garch_fit_dist_mod,   only: fit_garch_dist_model, garch_dist_variance_path, &
                                    garch_dist_fit_result_t, model_param_count, dist_param_count
    use msgarch_mod,          only: msgarch_result_t, fit_msgarch, &
                                    set_msgarch_analytic_grad, set_msgarch_gauss_noise
    use msnagarch_mod,        only: msnagarch_result_t, fit_msnagarch, &
                                    set_msnagarch_analytic_grad, set_msnagarch_gauss_noise
    implicit none

    integer,  parameter :: nmod = 4
    integer,  parameter :: max_iter = 500
    real(dp), parameter :: gtol = 1.0e-6_dp
    real(dp), parameter :: trading_days = 252.0_dp
    character(len=9), parameter :: model_names(nmod) = &
        [character(len=9) :: "GARCH", "NAGARCH", "MSGARCH", "MSNAGARCH"]
    integer, parameter :: ms_k_base(2) = [8, 10]   ! MSGARCH, MSNAGARCH param counts without dof

    character(len=256) :: prices_file
    integer,  allocatable :: dates(:)
    character(len=32), allocatable :: col_names(:)
    real(dp), allocatable :: prices(:,:), ret(:), h(:)

    type(garch_dist_fit_result_t) :: garch_res, nagarch_res
    type(msgarch_result_t)        :: ms_res
    type(msnagarch_result_t)      :: msn_res

    real(dp) :: logL(nmod), aic(nmod), bic(nmod), vol_ann(nmod), t_fit(nmod)
    real(dp) :: t0, t1
    integer  :: k(nmod), niters(nmod), aic_rank(nmod), bic_rank(nmod)
    integer  :: nprices, ncols, nobs, icol, nargs, id, i
    logical  :: conv(nmod), analytic_grad, gauss
    character(len=32) :: arg, dist_arg
    character(len=16) :: dist_names(4)
    integer  :: ndist

    call print_program_header("xcompare_ms_garch_models.f90")
    prices_file      = "spy_efa_eem_tlt_lqd.csv" ! "vix_spy.csv"
    analytic_grad = .true.
    dist_arg      = "NORMAL,T"

    nargs = command_argument_count()
    i = 1
    do while (i <= nargs)
        call get_command_argument(i, arg)
        if (trim(arg) == "-numerical") then
            analytic_grad = .false.
        else if (trim(arg) == "-dist") then
            i = i + 1
            if (i <= nargs) call get_command_argument(i, dist_arg)
        else if (index(trim(arg), ".csv") > 0 .or. index(trim(arg), ".CSV") > 0) then
            prices_file = arg
        else if (i == 1) then
            prices_file = arg
        end if
        i = i + 1
    end do
    call uppercase_str(dist_arg)
    call parse_csv(trim(dist_arg), dist_names, ndist)

    call set_msgarch_analytic_grad(analytic_grad)
    call set_msnagarch_analytic_grad(analytic_grad)
    write(*, '(A,L1)') "Gradient mode (analytic): ", analytic_grad

    call read_price_csv(trim(prices_file), dates, col_names, prices)
    nprices = size(prices, 1)
    ncols   = size(prices, 2)
    nobs    = nprices - 1
    allocate(ret(nobs), h(nobs))
    call print_price_sample_info(trim(prices_file), dates, ncols, nobs)

    do icol = 1, ncols
        ret = log(prices(2:nprices, icol) / prices(1:nprices-1, icol))
        ret = ret - mean(ret)

        do id = 1, ndist
            gauss = (trim(dist_names(id)) == "NORMAL")
            call set_msgarch_gauss_noise(gauss)
            call set_msnagarch_gauss_noise(gauss)

            ! Parameter counts: GARCH and NAGARCH via garch_fit_dist_mod;
            ! MS models drop dof for Gaussian.
            k(1) = model_param_count("GARCH")   + dist_param_count(trim(dist_names(id)))
            k(2) = model_param_count("NAGARCH")  + dist_param_count(trim(dist_names(id)))
            k(3) = ms_k_base(1) + merge(0, 1, gauss)
            k(4) = ms_k_base(2) + merge(0, 1, gauss)

            ! --- GARCH ---
            call cpu_time(t0)
            call fit_garch_dist_model("GARCH", trim(dist_names(id)), ret, max_iter, gtol, garch_res)
            call cpu_time(t1)
            t_fit(1)  = t1 - t0
            logL(1)   = -garch_res%nll * real(nobs, dp)
            niters(1) = garch_res%niter
            conv(1)   = garch_res%converged
            call garch_dist_variance_path("GARCH", ret, garch_res%params, h)
            vol_ann(1) = unc_vol(h, nobs)

            ! --- NAGARCH ---
            call cpu_time(t0)
            call fit_garch_dist_model("NAGARCH", trim(dist_names(id)), ret, max_iter, gtol, nagarch_res)
            call cpu_time(t1)
            t_fit(2)  = t1 - t0
            logL(2)   = -nagarch_res%nll * real(nobs, dp)
            niters(2) = nagarch_res%niter
            conv(2)   = nagarch_res%converged
            call garch_dist_variance_path("NAGARCH", ret, nagarch_res%params, h)
            vol_ann(2) = unc_vol(h, nobs)

            ! --- MSGARCH ---
            call cpu_time(t0)
            call fit_msgarch(ret, nobs, max_iter, gtol, ms_res, h)
            call cpu_time(t1)
            t_fit(3)  = t1 - t0
            logL(3)   = ms_res%loglik
            niters(3) = ms_res%niter
            conv(3)   = ms_res%converged
            vol_ann(3) = unc_vol(h, nobs)

            ! --- MSNAGARCH ---
            call cpu_time(t0)
            call fit_msnagarch(ret, nobs, max_iter, gtol, msn_res, h)
            call cpu_time(t1)
            t_fit(4)  = t1 - t0
            logL(4)   = msn_res%loglik
            niters(4) = msn_res%niter
            conv(4)   = msn_res%converged
            vol_ann(4) = unc_vol(h, nobs)

            aic = 2.0_dp*real(k, dp) - 2.0_dp*logL
            bic = log(real(nobs, dp))*real(k, dp) - 2.0_dp*logL
            call rank_asc(aic, aic_rank)
            call rank_asc(bic, bic_rank)

            call print_summary(col_names(icol), nobs, trim(dist_names(id)), &
                               logL, k, aic, bic, aic_rank, bic_rank, vol_ann, t_fit, niters, conv)
            call print_params(garch_res, nagarch_res, ms_res, msn_res, gauss)
            write(*, '(A)') ""
        end do
    end do

contains

    real(dp) function unc_vol(h, n)
        real(dp), intent(in) :: h(:)
        integer,  intent(in) :: n
        unc_vol = sqrt(trading_days * sum(h) / real(n, dp)) * 100.0_dp
    end function unc_vol

    subroutine rank_asc(x, r)
        real(dp), intent(in)  :: x(:)
        integer,  intent(out) :: r(:)
        integer :: j
        do j = 1, size(x)
            r(j) = 1 + count(x < x(j))
        end do
    end subroutine rank_asc

    subroutine uppercase_str(s)
        character(len=*), intent(inout) :: s
        integer :: j, c
        do j = 1, len(s)
            c = iachar(s(j:j))
            if (c >= iachar('a') .and. c <= iachar('z')) &
                s(j:j) = achar(c - 32)
        end do
    end subroutine uppercase_str

    subroutine parse_csv(s, names, n)
        character(len=*), intent(in)  :: s
        character(len=16), intent(out) :: names(:)
        integer,           intent(out) :: n
        integer :: pos, comma
        n   = 0
        pos = 1
        do while (pos <= len_trim(s) .and. n < size(names))
            comma = index(s(pos:), ",")
            if (comma == 0) then
                n = n + 1
                names(n) = adjustl(s(pos:))
                exit
            else
                n = n + 1
                names(n) = adjustl(s(pos:pos+comma-2))
                pos = pos + comma
            end if
        end do
    end subroutine parse_csv

    subroutine print_summary(asset, n, dist, logL_in, k_in, aic, bic, &
                             aic_rank, bic_rank, vol_ann, t_fit, niters, conv)
        character(len=*), intent(in) :: asset, dist
        integer,  intent(in) :: n, k_in(:), aic_rank(:), bic_rank(:), niters(:)
        real(dp), intent(in) :: logL_in(:), aic(:), bic(:), vol_ann(:), t_fit(:)
        logical,  intent(in) :: conv(:)
        integer :: j
        character(len=*), parameter :: hdr = &
            "------------------------------------------------------------------------------------" // &
            "------------------------------"

        write(*, '(A,A,A,I0,A,A,A)') "Asset: ", trim(asset), "   n = ", n, &
            "   distribution: ", trim(dist)
        write(*, '(A)') hdr
        write(*, '(A9,1X,A3,1X,A12,1X,A12,1X,A12,1X,A8,1X,A8,1X,A9,1X,A7,1X,A5,1X,A4)') &
            "Model", "k", "logL", "AIC", "BIC", "AIC_rank", "BIC_rank", "vol_ann%", "cpu(s)", "iter", "conv"
        write(*, '(A)') hdr
        do j = 1, nmod
            write(*, '(A9,1X,I3,1X,F12.3,1X,F12.3,1X,F12.3,1X,I8,1X,I8,1X,F9.3,1X,F7.2,1X,I5,1X,L4)') &
                model_names(j), k_in(j), logL_in(j), aic(j), bic(j), &
                aic_rank(j), bic_rank(j), vol_ann(j), t_fit(j), niters(j), conv(j)
        end do
        write(*, '(A)') hdr
    end subroutine print_summary

    subroutine print_params(gr, nr, ms, msn, gauss)
        type(garch_dist_fit_result_t), intent(in) :: gr, nr
        type(msgarch_result_t),        intent(in) :: ms
        type(msnagarch_result_t),      intent(in) :: msn
        logical,                       intent(in) :: gauss

        write(*, '(A)') ""
        write(*, '(A)') "GARCH parameters:"
        if (gauss) then
            write(*, '(2X,A,ES12.4,2X,A,F8.5,2X,A,F8.5,2X,A,F8.5)') &
                "omega=", gr%params%omega, "alpha=", gr%params%alpha, &
                "beta=", gr%params%beta, "persist=", gr%persist
        else
            write(*, '(2X,A,ES12.4,2X,A,F8.5,2X,A,F8.5,2X,A,F8.5,2X,A,F6.2)') &
                "omega=", gr%params%omega, "alpha=", gr%params%alpha, &
                "beta=", gr%params%beta, "persist=", gr%persist, "nu=", gr%shape
        end if
        write(*, '(A)') "NAGARCH parameters:"
        if (gauss) then
            write(*, '(2X,A,ES12.4,2X,A,F8.5,2X,A,F8.5,2X,A,F8.5,2X,A,F8.5)') &
                "omega=", nr%params%omega, "alpha=", nr%params%alpha, &
                "beta=", nr%params%beta, "theta=", nr%params%theta, "persist=", nr%persist
        else
            write(*, '(2X,A,ES12.4,2X,A,F8.5,2X,A,F8.5,2X,A,F8.5,2X,A,F8.5,2X,A,F6.2)') &
                "omega=", nr%params%omega, "alpha=", nr%params%alpha, &
                "beta=", nr%params%beta, "theta=", nr%params%theta, &
                "persist=", nr%persist, "nu=", nr%shape
        end if
        write(*, '(A)') "MSGARCH parameters:"
        if (gauss) then
            write(*, '(2X,A,F7.4,2X,A,F7.4,2X,A,F7.4,2X,A,F7.4,2X,A,F7.4,2X,A,F7.4)') &
                "persist1=", ms%persist(1), "persist2=", ms%persist(2), &
                "alpha1=", ms%params%alpha(1), "alpha2=", ms%params%alpha(2), &
                "p11=", ms%params%p11, "p22=", ms%params%p22
        else
            write(*, '(2X,A,F7.4,2X,A,F7.4,2X,A,F7.4,2X,A,F7.4,2X,A,F7.4,2X,A,F7.4,2X,A,F6.2)') &
                "persist1=", ms%persist(1), "persist2=", ms%persist(2), &
                "alpha1=", ms%params%alpha(1), "alpha2=", ms%params%alpha(2), &
                "p11=", ms%params%p11, "p22=", ms%params%p22, "dof=", ms%params%dof
        end if
        write(*, '(2X,A,ES11.3,2X,A,ES11.3)') &
            "omega1=", ms%params%omega(1), "omega2=", ms%params%omega(2)
        write(*, '(A)') "MSNAGARCH parameters:"
        if (gauss) then
            write(*, '(2X,A,F7.4,2X,A,F7.4,2X,A,F7.4,2X,A,F7.4,2X,A,F7.4,2X,A,F7.4)') &
                "persist1=", msn%persist(1), "persist2=", msn%persist(2), &
                "theta1=", msn%params%theta(1), "theta2=", msn%params%theta(2), &
                "p11=", msn%params%p11, "p22=", msn%params%p22
        else
            write(*, '(2X,A,F7.4,2X,A,F7.4,2X,A,F7.4,2X,A,F7.4,2X,A,F7.4,2X,A,F7.4,2X,A,F6.2)') &
                "persist1=", msn%persist(1), "persist2=", msn%persist(2), &
                "theta1=", msn%params%theta(1), "theta2=", msn%params%theta(2), &
                "p11=", msn%params%p11, "p22=", msn%params%p22, "dof=", msn%params%dof
        end if
        write(*, '(2X,A,ES11.3,2X,A,ES11.3)') &
            "omega1=", msn%params%omega(1), "omega2=", msn%params%omega(2)
    end subroutine print_params

end program xcompare_ms_garch_models

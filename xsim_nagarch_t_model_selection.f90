! Monte Carlo study of AIC/BIC model selection when the DGP is NAGARCH-t.
!
! For each replication:
!   1. Simulate nobs returns from NAGARCH with standardised Student-t innovations.
!   2. Fit all four models (GARCH, NAGARCH, MSGARCH, MSNAGARCH) twice: once
!      assuming Normal noise and once assuming Student-t noise.
!   3. Compute AIC and BIC for each model×distribution combination and record
!      which model wins each criterion.
!
! Usage: xsim_nagarch_t_model_selection [nrep] [nobs]
!   nrep : Monte Carlo replications (default 100)
!   nobs : observations per replication (default 2000)

program xsim_nagarch_t_model_selection
    use kind_mod,           only: dp
    use date_mod,           only: print_program_header
    use stats_mod,          only: mean
    use random_mod,         only: random_t_std
    use garch_fit_dist_mod, only: fit_garch_dist_model, garch_dist_fit_result_t, &
                                  model_param_count, dist_param_count
    use msgarch_mod,        only: msgarch_result_t, fit_msgarch, &
                                  set_msgarch_analytic_grad, set_msgarch_gauss_noise
    use msnagarch_mod,      only: msnagarch_result_t, fit_msnagarch, &
                                  set_msnagarch_analytic_grad, set_msnagarch_gauss_noise
    implicit none

    integer,  parameter :: nmod     = 4
    integer,  parameter :: ndist    = 2
    integer,  parameter :: max_iter = 500
    integer,  parameter :: nburn    = 500
    real(dp), parameter :: gtol     = 1.0e-6_dp
    ! DGP: NAGARCH-t
    real(dp), parameter :: dgp_omega = 2.5e-6_dp
    real(dp), parameter :: dgp_alpha = 0.08_dp
    real(dp), parameter :: dgp_beta  = 0.80_dp
    real(dp), parameter :: dgp_theta = 1.0_dp
    real(dp), parameter :: dgp_nu    = 7.0_dp
    integer,  parameter :: ms_k_base(2) = [8, 10]

    character(len=9), parameter :: model_names(nmod) = &
        [character(len=9) :: "GARCH", "NAGARCH", "MSGARCH", "MSNAGARCH"]
    character(len=6), parameter :: dist_names(ndist) = &
        [character(len=6) :: "NORMAL", "T"]

    real(dp), allocatable :: ret(:), h_ms(:)
    type(garch_dist_fit_result_t) :: garch_res, nagarch_res
    type(msgarch_result_t)        :: ms_res
    type(msnagarch_result_t)      :: msn_res

    integer  :: aic_wins(nmod, ndist), bic_wins(nmod, ndist), nconv(nmod, ndist)
    real(dp) :: sum_logL(nmod, ndist)
    real(dp) :: logL_rep(nmod), aic_rep(nmod), bic_rep(nmod)
    integer  :: k_arr(nmod)
    logical  :: gauss
    integer  :: nrep, nobs, irep, id, aic_winner, bic_winner, ios
    character(len=32) :: arg

    call print_program_header("xsim_nagarch_t_model_selection.f90")
    nrep = 100
    nobs = 2000

    if (command_argument_count() >= 1) then
        call get_command_argument(1, arg)
        read(arg, *, iostat=ios) nrep
        if (ios /= 0 .or. nrep < 1) error stop "invalid nrep"
    end if
    if (command_argument_count() >= 2) then
        call get_command_argument(2, arg)
        read(arg, *, iostat=ios) nobs
        if (ios /= 0 .or. nobs < 100) error stop "invalid nobs"
    end if

    allocate(ret(nobs), h_ms(nobs))
    call set_msgarch_analytic_grad(.true.)
    call set_msnagarch_analytic_grad(.true.)

    write(*, '(A)') "DGP: NAGARCH with Student-t innovations"
    write(*, '(2X,A,ES10.3,2X,A,F5.3,2X,A,F5.3,2X,A,F5.3,2X,A,F5.2)') &
        "omega=", dgp_omega, "alpha=", dgp_alpha, "beta=", dgp_beta, &
        "theta=", dgp_theta, "nu=", dgp_nu
    write(*, '(2X,A,F6.4)') "persistence=", dgp_alpha*(1.0_dp + dgp_theta**2) + dgp_beta
    write(*, '(A,I0,A,I0)') "Replications: ", nrep, "   nobs per rep: ", nobs
    write(*, '(A)') ""

    aic_wins = 0
    bic_wins = 0
    nconv    = 0
    sum_logL = 0.0_dp

    do irep = 1, nrep
        call seed_rng(irep)
        call simulate_nagarch_t(nobs, ret)
        ret = ret - mean(ret)

        do id = 1, ndist
            gauss = (trim(adjustl(dist_names(id))) == "NORMAL")
            call set_msgarch_gauss_noise(gauss)
            call set_msnagarch_gauss_noise(gauss)

            k_arr(1) = model_param_count("GARCH")   + dist_param_count(trim(adjustl(dist_names(id))))
            k_arr(2) = model_param_count("NAGARCH")  + dist_param_count(trim(adjustl(dist_names(id))))
            k_arr(3) = ms_k_base(1) + merge(0, 1, gauss)
            k_arr(4) = ms_k_base(2) + merge(0, 1, gauss)

            call fit_garch_dist_model("GARCH",   trim(adjustl(dist_names(id))), ret, max_iter, gtol, garch_res)
            call fit_garch_dist_model("NAGARCH",  trim(adjustl(dist_names(id))), ret, max_iter, gtol, nagarch_res)
            call fit_msgarch(ret, nobs, max_iter, gtol, ms_res, h_ms)
            call fit_msnagarch(ret, nobs, max_iter, gtol, msn_res, h_ms)

            logL_rep(1) = -garch_res%nll   * real(nobs, dp)
            logL_rep(2) = -nagarch_res%nll  * real(nobs, dp)
            logL_rep(3) = ms_res%loglik
            logL_rep(4) = msn_res%loglik

            if (garch_res%converged)   nconv(1, id) = nconv(1, id) + 1
            if (nagarch_res%converged) nconv(2, id) = nconv(2, id) + 1
            if (ms_res%converged)      nconv(3, id) = nconv(3, id) + 1
            if (msn_res%converged)     nconv(4, id) = nconv(4, id) + 1

            sum_logL(:, id) = sum_logL(:, id) + logL_rep

            aic_rep = 2.0_dp*real(k_arr, dp) - 2.0_dp*logL_rep
            bic_rep = log(real(nobs, dp))*real(k_arr, dp) - 2.0_dp*logL_rep

            aic_winner = minloc(aic_rep, 1)
            bic_winner = minloc(bic_rep, 1)
            aic_wins(aic_winner, id) = aic_wins(aic_winner, id) + 1
            bic_wins(bic_winner, id) = bic_wins(bic_winner, id) + 1
        end do

        if (mod(irep, max(1, nrep/10)) == 0) &
            write(*, '(A,I0,A,I0)') "  rep ", irep, " / ", nrep
    end do

    call print_results(nrep, nobs)

contains

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

    subroutine simulate_nagarch_t(n, y)
        integer,  intent(in)  :: n
        real(dp), intent(out) :: y(:)
        real(dp) :: h, sqrth, eps, persist
        integer  :: t
        persist = dgp_alpha*(1.0_dp + dgp_theta**2) + dgp_beta
        h = max(dgp_omega / max(1.0_dp - persist, 1.0e-8_dp), 1.0e-12_dp)
        do t = 1, nburn
            sqrth = sqrt(h)
            eps   = random_t_std(dgp_nu)
            h = max(dgp_omega + dgp_alpha*(sqrth*eps - dgp_theta*sqrth)**2 + dgp_beta*h, 1.0e-12_dp)
        end do
        do t = 1, n
            sqrth = sqrt(h)
            eps   = random_t_std(dgp_nu)
            y(t)  = sqrth * eps
            h = max(dgp_omega + dgp_alpha*(y(t) - dgp_theta*sqrth)**2 + dgp_beta*h, 1.0e-12_dp)
        end do
    end subroutine simulate_nagarch_t

    subroutine print_results(nrep, nobs)
        integer, intent(in) :: nrep, nobs
        integer  :: j, id, k_loc(nmod)
        logical  :: gauss_loc
        real(dp) :: fn
        character(len=*), parameter :: hdr = repeat("-", 76)

        fn = real(nrep, dp)
        write(*, '(/,A,I0,A,I0)') &
            "Model selection frequency: DGP=NAGARCH-t  nrep=", nrep, "  nobs=", nobs

        do id = 1, ndist
            gauss_loc = (trim(adjustl(dist_names(id))) == "NORMAL")
            k_loc(1) = model_param_count("GARCH")   + dist_param_count(trim(adjustl(dist_names(id))))
            k_loc(2) = model_param_count("NAGARCH")  + dist_param_count(trim(adjustl(dist_names(id))))
            k_loc(3) = ms_k_base(1) + merge(0, 1, gauss_loc)
            k_loc(4) = ms_k_base(2) + merge(0, 1, gauss_loc)

            write(*, '(A)') ""
            write(*, '(A,A)') "Fitted distribution: ", trim(adjustl(dist_names(id)))
            write(*, '(A)') hdr
            write(*, '(A9,1X,A3,1X,A8,1X,A6,1X,A8,1X,A6,1X,A6,1X,A12)') &
                "Model", "k", "AIC_wins", "AIC_%", "BIC_wins", "BIC_%", "conv_%", "mean_logL"
            write(*, '(A)') hdr
            do j = 1, nmod
                write(*, '(A9,1X,I3,1X,I8,1X,F6.1,1X,I8,1X,F6.1,1X,F6.1,1X,F12.2)') &
                    model_names(j), k_loc(j), &
                    aic_wins(j,id),  100.0_dp*aic_wins(j,id)/fn, &
                    bic_wins(j,id),  100.0_dp*bic_wins(j,id)/fn, &
                    100.0_dp*nconv(j,id)/fn, &
                    sum_logL(j,id)/fn
            end do
            write(*, '(A)') hdr
        end do
    end subroutine print_results

end program xsim_nagarch_t_model_selection

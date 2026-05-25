! Simulate from GAS(1,1)-Normal and fit both GAS(1,1)-Normal and GARCH(1,1)-Normal.
! Since data are GAS-generated, GAS should achieve higher log-likelihoods.
! Branched from xgas_scaling.f90.

program xgas_garch_scaling
    use kind_mod,       only: dp
    use gas_mod,        only: gas_simulate, gas_set_data, gas_set_types, gas_np, gas_obj, &
                              gas_sym_inv_transform, gas_transform, &
                              proc_gas, dist_normal
    use garch_flex_mod, only: flex_set_data, flex_set_types, flex_np, flex_obj, proc_garch
    use garch_mod,   only: garch_inv_transform, garch_transform
    use bfgs_mod,    only: bfgs_minimize
    implicit none

    ! ---- tuneable parameters ----
    integer,  parameter :: n_sizes  = 6
    integer,  parameter :: n_start  = 250
    integer,  parameter :: n_mult   = 4
    integer,  parameter :: seed_val = 42
    integer,  parameter :: max_iter = 500
    real(dp), parameter :: gtol     = 1.0e-7_dp

    ! ---- true GAS(1,1)-Normal parameters ----
    real(dp), parameter :: true_omega = 0.0_dp
    real(dp), parameter :: true_alpha = 0.10_dp
    real(dp), parameter :: true_beta  = 0.95_dp

    real(dp), allocatable :: y(:), p(:), p0(:)
    real(dp) :: omega, alpha, gamma_v, beta, f_opt
    logical  :: converged
    integer  :: np, n_iter, i, n

    integer  :: sizes(n_sizes)
    real(dp) :: gas_omega(n_sizes),   gas_alpha(n_sizes),   gas_beta(n_sizes),   gas_logl(n_sizes)
    real(dp) :: garch_omega(n_sizes), garch_alpha(n_sizes), garch_beta(n_sizes), garch_logl(n_sizes)

    ! GAS row: n, model, omega, alpha, beta, logL
    character(len=*), parameter :: gas_fmt   = "(I9, A8, ES11.3, F8.4, F9.4, F12.1)"
    ! GARCH row: (blank n), model, omega, alpha, beta, logL, GAS-GARCH delta
    character(len=*), parameter :: garch_fmt = "(9X, A8, ES11.3, F8.4, F9.4, F12.1, F9.1)"

    n = n_start
    do i = 1, n_sizes
        sizes(i) = n
        n = n * n_mult
    end do

    do i = 1, n_sizes
        n = sizes(i)
        allocate(y(n))
        call gas_simulate(true_omega, true_alpha, 0.0_dp, true_beta, &
                          dist_normal, 0.0_dp, n, seed_val, y)

        ! ── fit GAS(1,1)-Normal ──────────────────────────────────────────────
        call gas_set_data(y, n)
        call gas_set_types(proc_gas, dist_normal)
        np = gas_np()
        allocate(p(np), p0(np))
        call gas_sym_inv_transform(0.0_dp, 0.05_dp, 0.90_dp, p0)
        p = p0
        call bfgs_minimize(gas_obj, p, np, max_iter, gtol, f_opt, n_iter, converged)
        call gas_transform(p, omega, alpha, gamma_v, beta)
        gas_omega(i) = omega
        gas_alpha(i) = alpha
        gas_beta(i)  = beta
        gas_logl(i)  = -n * f_opt
        deallocate(p, p0)

        ! ── fit GARCH(1,1)-Normal ────────────────────────────────────────────
        call flex_set_data(y, n)
        call flex_set_types(proc_garch, dist_normal)
        np = flex_np()
        allocate(p(np), p0(np))
        call garch_inv_transform(1.0e-5_dp, 0.08_dp, 0.88_dp, p0)
        p = p0
        call bfgs_minimize(flex_obj, p, np, max_iter, gtol, f_opt, n_iter, converged)
        call garch_transform(p(1:3), omega, alpha, beta)
        garch_omega(i) = omega
        garch_alpha(i) = alpha
        garch_beta(i)  = beta
        garch_logl(i)  = -n * f_opt
        deallocate(p, p0)

        deallocate(y)
    end do

    print '(A)', ""
    print '(A)', " GAS(1,1)-Normal data: GAS vs GARCH fit comparison"
    print '(A,I0,A,I0,A)', " Sizes: ", n_start, " x ", n_mult, "^(0,1,2,...)"
    print '(A)', " True: omega=0.00  alpha=0.10  beta=0.95"
    print '(A)', ""
    print '(A9, A8, A11, A8, A9, A12, A9)', &
        "n       ", "Model   ", "omega      ", "alpha   ", "beta     ", "logL        ", "GAS-GARCH"
    print '(A)', repeat("-", 66)
    do i = 1, n_sizes
        write(*, gas_fmt)   sizes(i), "GAS     ", gas_omega(i),   gas_alpha(i),   gas_beta(i),   gas_logl(i)
        write(*, garch_fmt)           "GARCH   ", garch_omega(i), garch_alpha(i), garch_beta(i), garch_logl(i), &
            gas_logl(i) - garch_logl(i)
        write(*, *)
    end do

end program xgas_garch_scaling

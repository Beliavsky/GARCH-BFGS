! Fit GAS(1,1)-Normal on successively larger samples from the same process.
! Sample sizes are n_start, n_start*n_mult, n_start*n_mult^2, ...
! Prints a timing and estimation summary at the end.

program xgas_scaling
    use kind_mod,    only: dp
    use gas_mod,     only: gas_simulate, gas_set_data, gas_set_types, gas_np, gas_obj, &
                           gas_sym_inv_transform, gas_transform, &
                           proc_gas, dist_normal
    use bfgs_mod, only: bfgs_minimize
    implicit none

    ! ---- tuneable parameters ----
    integer,  parameter :: n_sizes  = 3
    integer,  parameter :: n_start  = 250
    integer,  parameter :: n_mult   = 4
    integer,  parameter :: seed_val = 42
    integer,  parameter :: max_iter = 500
    real(dp), parameter :: gtol     = 1.0e-7_dp

    ! ---- true GAS(1,1)-Normal parameters ----
    real(dp), parameter :: true_omega = 0.0_dp    ! log-variance intercept
    real(dp), parameter :: true_alpha = 0.10_dp   ! score coefficient
    real(dp), parameter :: true_beta  = 0.95_dp   ! persistence

    real(dp), allocatable :: y(:), p(:), p0(:)
    real(dp) :: omega, alpha, gamma_v, beta, f_opt
    logical  :: converged
    integer  :: np, n_iter, i, n

    integer  :: sizes(n_sizes)
    integer  :: ck0, ck1, ck_rate
    real(dp) :: t_sim(n_sizes), t_opt(n_sizes)
    real(dp) :: est_omega(n_sizes), est_alpha(n_sizes), est_beta(n_sizes)
    integer  :: iters(n_sizes)

    call system_clock(count_rate=ck_rate)
    call gas_set_types(proc_gas, dist_normal)
    np = gas_np()

    n = n_start
    do i = 1, n_sizes
        sizes(i) = n
        n = n * n_mult
    end do

    do i = 1, n_sizes
        n = sizes(i)
        allocate(y(n), p(np), p0(np))

        call system_clock(count=ck0)
        call gas_simulate(true_omega, true_alpha, 0.0_dp, true_beta, &
                          dist_normal, 0.0_dp, n, seed_val, y)
        call gas_set_data(y, n)
        call system_clock(count=ck1)
        t_sim(i) = real(ck1 - ck0, dp) / real(ck_rate, dp)

        call gas_sym_inv_transform(0.0_dp, 0.05_dp, 0.90_dp, p0)
        p = p0
        call system_clock(count=ck0)
        call bfgs_minimize(gas_obj, p, np, max_iter, gtol, f_opt, n_iter, converged)
        call system_clock(count=ck1)
        t_opt(i) = real(ck1 - ck0, dp) / real(ck_rate, dp)

        call gas_transform(p, omega, alpha, gamma_v, beta)
        est_omega(i) = omega
        est_alpha(i) = alpha
        est_beta(i)  = beta
        iters(i)     = n_iter

        deallocate(y, p, p0)
    end do

    print '(A)', ""
    print '(A)', " GAS(1,1)-Normal scaling study"
    print '(A,I0,A,I0,A)', " Sizes: ", n_start, " x ", n_mult, "^(0,1,2,...)"
    print '(A)', " True: omega=0.00  alpha=0.10  beta=0.95"
    print '(A)', ""
    print '(A,I0)', " Max iterations: ", max_iter
    print '(A)', ""
    print '(A)', "        n    omega    alpha     beta iters      t_sim      t_opt      t_tot"
    print '(A)', "  -------  -------  -------  ------- -----  ---------  ---------  ---------"
    do i = 1, n_sizes
        print '(I9, 3(F9.4), I6, 3(F11.6))', &
            sizes(i), est_omega(i), est_alpha(i), est_beta(i), iters(i), &
            t_sim(i), t_opt(i), t_sim(i) + t_opt(i)
    end do
    print '(A)', ""

end program xgas_scaling

! Fit GARCH(1,1) on successively larger samples from the same process.
! Sample sizes are n_start, n_start*n_mult, n_start*n_mult^2, ...
! Prints a timing and estimation summary at the end.

program garch_scaling
    use kind_mod,     only: dp
    use garch_module, only: garch_simulate, garch_set_data, garch_inv_transform, &
                            garch_obj, garch_transform
    use bfgs_module,  only: bfgs_minimize
    implicit none

    ! ---- tuneable parameters ----
    integer,  parameter :: n_sizes  = 6
    integer,  parameter :: n_start  = 250
    integer,  parameter :: n_mult   = 4
    integer,  parameter :: seed_val = 42
    integer,  parameter :: max_iter = 500
    real(dp), parameter :: gtol     = 1.0e-7_dp

    ! ---- true GARCH(1,1) parameters ----
    real(dp), parameter :: true_omega = 0.05_dp
    real(dp), parameter :: true_alpha = 0.10_dp
    real(dp), parameter :: true_beta  = 0.85_dp

    real(dp), allocatable :: y(:)
    real(dp) :: p(3), omega, alpha, beta, f_opt
    logical  :: converged
    integer  :: n_iter, i, n

    integer  :: sizes(n_sizes)
    integer  :: ck0, ck1, ck_rate
    real(dp) :: t_sim(n_sizes), t_opt(n_sizes)
    real(dp) :: est_omega(n_sizes), est_alpha(n_sizes), est_beta(n_sizes)
    integer  :: iters(n_sizes)

    call system_clock(count_rate=ck_rate)

    ! ---- build size sequence ----
    n = n_start
    do i = 1, n_sizes
        sizes(i) = n
        n = n * n_mult
    end do

    ! ---- simulate, set data, and fit for each size ----
    do i = 1, n_sizes
        n = sizes(i)
        allocate(y(n))

        call system_clock(count=ck0)
        call garch_simulate(true_omega, true_alpha, true_beta, n, seed_val, y)
        call garch_set_data(y, n)
        call system_clock(count=ck1)
        t_sim(i) = real(ck1 - ck0, dp) / real(ck_rate, dp)

        call garch_inv_transform(0.05_dp, 0.10_dp, 0.80_dp, p)

        call system_clock(count=ck0)
        call bfgs_minimize(garch_obj, p, 3, max_iter, gtol, f_opt, n_iter, converged)
        call system_clock(count=ck1)
        t_opt(i) = real(ck1 - ck0, dp) / real(ck_rate, dp)

        call garch_transform(p, omega, alpha, beta)
        est_omega(i) = omega
        est_alpha(i) = alpha
        est_beta(i)  = beta
        iters(i)     = n_iter

        deallocate(y)
    end do

    ! ---- summary ----
    print '(A)', ""
    print '(A)', " GARCH(1,1) scaling study"
    print '(A,I0,A,I0,A)', " Sizes: ", n_start, " x ", n_mult, "^(0,1,2,...)"
    print '(A)', " True: omega=0.05  alpha=0.10  beta=0.85"
    print '(A)', ""
    print '(A,I0)', " Max iterations: ", max_iter
    print '(A)', ""
    ! column widths: I9  F9.4  F9.4  F9.4  I6  F11.6  F11.6  F11.6
    print '(A)', "        n    omega    alpha     beta iters      t_sim      t_opt      t_tot"
    print '(A)', "  -------  -------  -------  ------- -----  ---------  ---------  ---------"
    do i = 1, n_sizes
        print '(I9, 3(F9.4), I6, 3(F11.6))', &
            sizes(i), est_omega(i), est_alpha(i), est_beta(i), iters(i), &
            t_sim(i), t_opt(i), t_sim(i) + t_opt(i)
    end do
    print '(A)', ""

end program garch_scaling

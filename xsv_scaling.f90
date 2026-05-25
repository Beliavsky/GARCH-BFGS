! Fit log-Normal SV(1) on successively larger samples from the same process.
! Sample sizes are n_start, n_start*n_mult, n_start*n_mult^2, ...
! Prints a timing and estimation summary at the end.

program xsv_scaling
    use kind_mod,    only: dp
    use sv_mod,      only: sv_simulate, sv_set_data, sv_set_types, sv_np, sv_obj, &
                           sv_sym_inv_transform, sv_transform, proc_sv, dist_normal
    use bfgs_mod, only: bfgs_minimize
    implicit none

    ! ---- tuneable parameters ----
    integer,  parameter :: n_sizes  = 6
    integer,  parameter :: n_start  = 250
    integer,  parameter :: n_mult   = 4
    integer,  parameter :: seed_val = 42
    integer,  parameter :: max_iter = 500
    real(dp), parameter :: gtol     = 1.0e-7_dp

    ! ---- true SV(1)-Normal parameters ----
    real(dp), parameter :: true_mu        = 0.0_dp    ! unconditional log-variance
    real(dp), parameter :: true_phi       = 0.97_dp   ! persistence
    real(dp), parameter :: true_sigma_eta = 0.15_dp   ! vol of log-volatility

    real(dp), allocatable :: y(:), p(:), p0(:)
    real(dp) :: mu, phi, sigma_eta, rho, nu, f_opt
    logical  :: converged
    integer  :: np, n_iter, i, n

    integer  :: sizes(n_sizes)
    integer  :: ck0, ck1, ck_rate
    real(dp) :: t_sim(n_sizes), t_opt(n_sizes)
    real(dp) :: est_mu(n_sizes), est_phi(n_sizes), est_se(n_sizes)
    integer  :: iters(n_sizes)

    call system_clock(count_rate=ck_rate)
    call sv_set_types(proc_sv, dist_normal)
    np = sv_np()

    n = n_start
    do i = 1, n_sizes
        sizes(i) = n
        n = n * n_mult
    end do

    do i = 1, n_sizes
        n = sizes(i)
        allocate(y(n), p(np), p0(np))

        call system_clock(count=ck0)
        call sv_simulate(true_mu, true_phi, true_sigma_eta, 0.0_dp, dist_normal, 0.0_dp, n, seed_val, y)
        call sv_set_data(y, n)
        call system_clock(count=ck1)
        t_sim(i) = real(ck1 - ck0, dp) / real(ck_rate, dp)

        call sv_sym_inv_transform(0.0_dp, 0.90_dp, 0.30_dp, p0)
        p = p0
        call system_clock(count=ck0)
        call bfgs_minimize(sv_obj, p, np, max_iter, gtol, f_opt, n_iter, converged)
        call system_clock(count=ck1)
        t_opt(i) = real(ck1 - ck0, dp) / real(ck_rate, dp)

        call sv_transform(p, mu, phi, sigma_eta, rho, nu)
        est_mu(i) = mu
        est_phi(i) = phi
        est_se(i)  = sigma_eta
        iters(i)   = n_iter

        deallocate(y, p, p0)
    end do

    print '(A)', ""
    print '(A)', " SV(1)-Normal QML scaling study"
    print '(A,I0,A,I0,A)', " Sizes: ", n_start, " x ", n_mult, "^(0,1,2,...)"
    print '(A)', " True: mu=0.00  phi=0.97  sigma_eta=0.15"
    print '(A)', ""
    print '(A,I0)', " Max iterations: ", max_iter
    print '(A)', ""
    print '(A)', "        n       mu      phi  sig_eta iters      t_sim      t_opt      t_tot"
    print '(A)', "  -------  -------  -------  ------- -----  ---------  ---------  ---------"
    do i = 1, n_sizes
        print '(I9, 3(F9.4), I6, 3(F11.6))', &
            sizes(i), est_mu(i), est_phi(i), est_se(i), iters(i), &
            t_sim(i), t_opt(i), t_sim(i) + t_opt(i)
    end do
    print '(A)', ""

end program xsv_scaling

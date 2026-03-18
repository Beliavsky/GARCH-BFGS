! Simulate from SV(1)-Normal with leverage and fit both SV(1) and SV(1)-lev.
! Since data are SV-lev-generated, SV-lev should achieve higher log-likelihoods.
! Branched from xsv_scaling.f90.

program xsv_lev_scaling
    use kind_mod,    only: dp
    use sv_mod,      only: sv_simulate, sv_set_data, sv_set_types, sv_np, sv_obj, &
                           sv_sym_inv_transform, sv_lev_inv_transform, sv_transform, &
                           proc_sv, proc_sv_lev, dist_normal
    use bfgs_module, only: bfgs_minimize
    implicit none

    ! ---- tuneable parameters ----
    integer,  parameter :: n_sizes  = 6
    integer,  parameter :: n_start  = 250
    integer,  parameter :: n_mult   = 4
    integer,  parameter :: seed_val = 42
    integer,  parameter :: max_iter = 500
    real(dp), parameter :: gtol     = 1.0e-7_dp

    ! ---- true SV(1)-lev-Normal parameters ----
    real(dp), parameter :: true_mu        =  0.0_dp   ! unconditional log-variance
    real(dp), parameter :: true_phi       =  0.97_dp  ! persistence
    real(dp), parameter :: true_sigma_eta =  0.15_dp  ! vol of log-volatility
    real(dp), parameter :: true_rho       = -0.6_dp   ! leverage correlation

    real(dp), allocatable :: y(:), p(:), p0(:)
    real(dp) :: mu, phi, sigma_eta, rho, nu, f_opt
    logical  :: converged
    integer  :: np, n_iter, i, n

    integer  :: sizes(n_sizes)
    real(dp) :: sv_mu(n_sizes),  sv_phi(n_sizes),  sv_se(n_sizes),  sv_logl(n_sizes)
    real(dp) :: lev_mu(n_sizes), lev_phi(n_sizes), lev_se(n_sizes), lev_rho(n_sizes), lev_logl(n_sizes)

    ! SV row:     n, model, mu, phi, sigma_eta, logL
    character(len=*), parameter :: sv_fmt  = "(I9, A9, F8.4, F8.4, F9.4, F12.1)"
    ! SV-lev row: (blank n), model, mu, phi, sigma_eta, rho, logL, SV-lev-SV delta
    character(len=*), parameter :: lev_fmt = "(9X, A9, F8.4, F8.4, F9.4, F8.4, F12.1, F9.1)"

    n = n_start
    do i = 1, n_sizes
        sizes(i) = n
        n = n * n_mult
    end do

    do i = 1, n_sizes
        n = sizes(i)
        allocate(y(n))
        call sv_simulate(true_mu, true_phi, true_sigma_eta, true_rho, dist_normal, 0.0_dp, n, seed_val, y)
        call sv_set_data(y, n)

        ! ── fit SV(1)-Normal ─────────────────────────────────────────────────
        call sv_set_types(proc_sv, dist_normal)
        np = sv_np()
        allocate(p(np), p0(np))
        call sv_sym_inv_transform(0.0_dp, 0.90_dp, 0.30_dp, p0)
        p = p0
        call bfgs_minimize(sv_obj, p, np, max_iter, gtol, f_opt, n_iter, converged)
        call sv_transform(p, mu, phi, sigma_eta, rho, nu)
        sv_mu(i)   = mu
        sv_phi(i)  = phi
        sv_se(i)   = sigma_eta
        sv_logl(i) = -n * f_opt
        deallocate(p, p0)

        ! ── fit SV(1)-lev-Normal ─────────────────────────────────────────────
        call sv_set_types(proc_sv_lev, dist_normal)
        np = sv_np()
        allocate(p(np), p0(np))
        call sv_lev_inv_transform(0.0_dp, 0.90_dp, 0.30_dp, -0.3_dp, p0)
        p = p0
        call bfgs_minimize(sv_obj, p, np, max_iter, gtol, f_opt, n_iter, converged)
        call sv_transform(p, mu, phi, sigma_eta, rho, nu)
        lev_mu(i)   = mu
        lev_phi(i)  = phi
        lev_se(i)   = sigma_eta
        lev_rho(i)  = rho
        lev_logl(i) = -n * f_opt
        deallocate(p, p0)

        deallocate(y)
    end do

    print '(A)', ""
    print '(A)', " SV(1)-lev-Normal data: SV vs SV-lev fit comparison"
    print '(A,I0,A,I0,A)', " Sizes: ", n_start, " x ", n_mult, "^(0,1,2,...)"
    print '(A)', " True: mu=0.00  phi=0.97  sigma_eta=0.15  rho=-0.60"
    print '(A)', ""
    print '(A9, A9, A8, A8, A9, A8, A12, A9)', &
        "n        ", "Model    ", "mu      ", "phi     ", "sig_eta  ", "rho     ", "logL        ", "lev-SV   "
    print '(A)', repeat("-", 72)
    do i = 1, n_sizes
        write(*, sv_fmt)  sizes(i), "SV       ", sv_mu(i),  sv_phi(i),  sv_se(i),              sv_logl(i)
        write(*, lev_fmt)           "SV-lev   ", lev_mu(i), lev_phi(i), lev_se(i), lev_rho(i), lev_logl(i), &
            lev_logl(i) - sv_logl(i)
        write(*, *)
    end do

end program xsv_lev_scaling

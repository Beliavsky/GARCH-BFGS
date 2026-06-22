! Standardised probability distributions: PDFs and iid fitting.
!
! All PDFs are unit-variance standardised forms (zero mean, unit variance),
! matching the innovation parameterisations used in the GARCH modules.
! For raw data y with mean mu and scale sigma, use x = (y-mu)/sigma.
!
! Distributions:
!   dist_normal   (1)  Normal                   no shape parameter
!   dist_t        (2)  Student-t, nu > 2         nu: degrees of freedom
!   dist_ged      (3)  GED, nu > 0               nu: shape (2=Normal, 1=Laplace)
!   dist_logistic (4)  Logistic                  no shape parameter
!   dist_laplace  (5)  Laplace                   no shape parameter
!   dist_sech     (6)  Hyperbolic secant         no shape parameter
!   dist_nig_sym  (7)  Symmetric NIG, alp > 0    alp: tail/shape parameter
!   dist_fs_skewt (8)  Fernandez-Steel skew-t    nu: df, xi: skew (>0)
!   dist_nig_gen  (9)  General NIG               alp: tail, rho=beta/alp: skew
!
! Additional PDFs (not wired into fit_dist):
!   pdf_nig_gen(x, alp, bet)  General NIG, alp > |bet|.
!     gam = sqrt(alp^2-bet^2),  delta = gam^3/alp^2,  mu0 = -bet*delta/gam.
!     Reduces to pdf_nig(x,alp) when bet=0.
!
!   pdf_vg_sym(x, nu)         Symmetric VG, nu > 0.
!     X = sqrt(G)*Z,  G ~ Gamma(1/nu, nu);  zero mean, unit variance.
!     f(x) proportional to K_{1/nu-1/2}(|x|*sqrt(2/nu)) * |x|^(1/nu-1/2).
!
!   pdf_vg_gen(x, nu, rho)    General (skewed) VG, nu > 0, rho in (-1,1).
!     rho = theta*sqrt(nu),  sigma^2 = 1-rho^2;  zero mean, unit variance.
!     X = theta*(G-1) + sigma*sqrt(G)*Z,  G ~ Gamma(1/nu, nu).
!     rho=0 recovers pdf_vg_sym(x, nu).
!
! Fitting routines:
!   fit_dist_std(x, n, dist_id, shape, loglik, converged)
!     Fit shape to standardised iid data x(1:n): mu=0 and sigma=1 fixed.
!     For no-shape distributions: NLL computed directly (converged=.true., shape=0).
!     For t, ged, nig: 1-parameter BFGS on the shape.
!     loglik = total (not per-obs) log-likelihood at solution.
!     Optional stats_start_in=.true. uses sample excess kurtosis for shape starts.
!     Optional robust_start_in=.true. uses quantile tail ratios for shape starts.
!
!   fit_dist(y, n, dist_id, mu, sigma, shape, loglik, converged)
!     Fit mu (= sample mean, fixed), sigma and shape to raw iid data y(1:n).
!     No-shape distributions: 1-parameter BFGS on sigma.
!     t, ged, nig: 2-parameter BFGS on (sigma, shape).
!     Optional stats_start_in=.true. uses sample excess kurtosis for shape starts.
!     Optional robust_start_in=.true. uses quantile tail ratios for shape starts.
!
! Parameter counts for fit_dist_std (useful for AIC/BIC of shape comparisons):
!   dist_npar_std = [0, 1, 1, 0, 0, 0, 1, 2]
!
! Score routines:
!   logpdf_std(x, dist_id, shape)       Standardised log-density.
!   score_z_std(x, dist_id, shape)      d log f(z;shape) / dz.
!   score_shape_std(x, dist_id, shape)  d log f(z;shape) / d shape.
! These are reusable in GARCH, GAS, SV, and iid distribution diagnostics.

module distributions_mod
    use kind_mod,       only: dp
    use math_const_mod, only: pi, two_pi, sqrt2, sqrt3, log2, half_log2, log_sqrt_2pi
    use bfgs_mod,    only: bfgs_minimize
    use special_mod,    only: bessel_k01, bessel_k_nu, digamma
    implicit none
    private

    ! ── Distribution ID constants ─────────────────────────────────────────────

    integer, parameter, public :: dist_normal   = 1
    integer, parameter, public :: dist_t        = 2
    integer, parameter, public :: dist_ged      = 3
    integer, parameter, public :: dist_logistic = 4
    integer, parameter, public :: dist_laplace  = 5
    integer, parameter, public :: dist_sech     = 6
    integer, parameter, public :: dist_nig_sym  = 7
    integer, parameter, public :: dist_nig      = dist_nig_sym ! compatibility alias for symmetric NIG
    integer, parameter, public :: dist_fs_skewt = 8
    integer, parameter, public :: dist_nig_gen  = 9
    integer, parameter, public :: dist_count    = 9

    character(len=8), dimension(dist_count), parameter, public :: dist_names = &
        ["normal  ", "t       ", "ged     ", "logistic", &
         "laplace ", "sech    ", "nig_sym ", "fs_skewt", "nig     "]

    ! Number of shape parameters estimated by fit_dist_std (0 or 1 per distribution).
    integer, dimension(dist_count), parameter, public :: dist_npar_std = &
        [0, 1, 1, 0, 0, 0, 1, 2, 2]

    ! ── Module-level state for BFGS callbacks ─────────────────────────────────

    real(dp), allocatable, save :: fd_y(:)     ! stored (demeaned) data
    integer,               save :: fd_n    = 0
    integer,               save :: fd_dist = 0
    integer,               save :: fd_mode = 0 ! 0 = sigma+shape fit, 1 = shape-only fit
    logical,               save :: fd_fixed_shape = .false.
    real(dp),              save :: fd_fixed_shape_val = 0.0_dp

    ! ── Public interface ──────────────────────────────────────────────────────

    public :: pdf_normal, pdf_t, pdf_ged, pdf_logistic, pdf_laplace, pdf_sech, pdf_nig
    public :: pdf_fs_skewt
    public :: pdf_nig_gen
    public :: pdf_vg_sym, pdf_vg_gen
    public :: logpdf_std, nll_std_val, score_z_std, score_shape_std
    public :: fit_dist_std, fit_dist
    public :: ged_lambda, dist_id_from_name, dist_fixed_shape_from_name
    public :: dist_warm_shape_start, dist_warm_shape_start_scaled, dist_warm_sigma_start

contains

    ! ── PDFs (all: zero mean, unit variance) ──────────────────────────────────

    pure elemental function pdf_normal(x) result(f)
        ! N(0,1): f(x) = (1/sqrt(2*pi)) * exp(-x^2/2)
        real(dp), intent(in) :: x
        real(dp) :: f
        f = exp(-0.5_dp*x**2 - log_sqrt_2pi)
    end function pdf_normal

    pure elemental function pdf_t(x, nu) result(f)
        ! Standardised t(nu), unit variance (nu > 2).
        ! f(x) = Γ((ν+1)/2) / (Γ(ν/2) * sqrt(π(ν-2))) * (1 + x²/(ν-2))^(-(ν+1)/2)
        real(dp), intent(in) :: x, nu
        real(dp) :: f
        f = exp(log_gamma(0.5_dp*(nu+1.0_dp)) - log_gamma(0.5_dp*nu) &
                - 0.5_dp*log(pi*(nu-2.0_dp)) &
                - 0.5_dp*(nu+1.0_dp)*log(1.0_dp + x**2/(nu-2.0_dp)))
    end function pdf_t

    pure elemental function pdf_ged(x, nu) result(f)
        ! Standardised GED(nu), unit variance.
        ! f(x) = ν/(2^(1+1/ν)*λ*Γ(1/ν)) * exp(-1/2*(|x|/λ)^ν),  λ = ged_lambda(ν)
        real(dp), intent(in) :: x, nu
        real(dp) :: f, lam
        lam = ged_lambda(nu)
        f   = exp(log(nu) - log(2.0_dp) &
                  - 1.5_dp*log_gamma(1.0_dp/nu) + 0.5_dp*log_gamma(3.0_dp/nu) &
                  - 0.5_dp*(abs(x)/lam)**nu)
    end function pdf_ged

    pure elemental function pdf_logistic(x) result(f)
        ! Standardised logistic, unit variance; c = pi/sqrt(3).
        ! f(x) = c * exp(-c|x|) / (1 + exp(-c|x|))^2
        real(dp), intent(in) :: x
        real(dp) :: f, c, au
        c  = pi / sqrt3
        au = c * abs(x)
        f  = c * exp(-au) / (1.0_dp + exp(-au))**2
    end function pdf_logistic

    pure elemental function pdf_laplace(x) result(f)
        ! Standardised Laplace, unit variance.
        ! f(x) = (1/sqrt(2)) * exp(-sqrt(2)*|x|)
        real(dp), intent(in) :: x
        real(dp) :: f
        f = exp(-sqrt2 * abs(x)) / sqrt2
    end function pdf_laplace

    pure elemental function pdf_sech(x) result(f)
        ! Standardised sech, unit variance.
        ! f(x) = (1/2) * sech(pi*x/2)
        real(dp), intent(in) :: x
        real(dp) :: f
        f = 0.5_dp / cosh(pi * x / 2.0_dp)
    end function pdf_sech

    pure elemental function pdf_nig(x, alp) result(f)
        ! Symmetric NIG, unit variance, alp > 0.
        ! f(x) = (alp^2/pi) * exp(alp^2) * K1(alp*sqrt(alp^2+x^2)) / sqrt(alp^2+x^2)
        real(dp), intent(in) :: x, alp
        real(dp) :: f, lk1, ratio, a2px2
        a2px2 = alp**2 + x**2
        call bessel_k01(alp * sqrt(a2px2), lk1, ratio)
        f = exp(2.0_dp*log(alp) - log(pi) + alp**2 + lk1 - 0.5_dp*log(a2px2))
    end function pdf_nig

    pure elemental function pdf_fs_skewt(x, nu, xi) result(f)
        ! Fernandez-Steel skewed Student-t, standardised to mean 0 and variance 1.
        ! xi=1 recovers the standardised symmetric t(nu).
        real(dp), intent(in) :: x, nu, xi
        real(dp) :: f, xip, xim, c, m1, raw_mean, raw_second, raw_sd, y, base_arg

        xip = max(xi, 1.0e-8_dp)
        xim = 1.0_dp / xip
        m1 = sqrt(nu - 2.0_dp) * exp(log_gamma(0.5_dp*(nu - 1.0_dp)) - &
             0.5_dp*log(pi) - log_gamma(0.5_dp*nu))
        raw_mean = m1 * (xip - xim)
        raw_second = (xip**3 + xim**3) / (xip + xim)
        raw_sd = sqrt(max(raw_second - raw_mean**2, 1.0e-12_dp))
        y = raw_mean + raw_sd*x
        c = 2.0_dp / (xip + xim)
        if (y >= 0.0_dp) then
            base_arg = y / xip
        else
            base_arg = y * xip
        end if
        f = raw_sd * c * pdf_t(base_arg, nu)
    end function pdf_fs_skewt

    pure elemental function pdf_nig_gen(x, alp, bet) result(f)
        ! General NIG, unit variance.  alp > |bet| >= 0.
        ! gam = sqrt(alp^2-bet^2),  delta = gam^3/alp^2,  mu0 = -bet*delta/gam.
        ! f(x) = (alp*delta/pi)*exp(delta*gam + bet*z)*K1(r)/sqrt(delta^2+z^2)
        !   where z = x - mu0,  r = alp*sqrt(delta^2+z^2).
        ! bet=0 recovers pdf_nig(x,alp) exactly.
        real(dp), intent(in) :: x, alp, bet   ! alp > |bet|
        real(dp) :: f, alp2, gam2, delta, mu0, dg, z, d2pz2, r, lk1, ratio
        alp2  = alp**2
        gam2  = alp2 - bet**2                 ! gamma^2
        delta = gam2 * sqrt(gam2) / alp2      ! gamma^3 / alpha^2
        mu0   = -bet * gam2 / alp2            ! zero-mean centering
        dg    = gam2 * gam2 / alp2            ! delta*gamma = gamma^4/alpha^2
        z     = x - mu0
        d2pz2 = delta**2 + z**2
        r     = alp * sqrt(d2pz2)
        call bessel_k01(r, lk1, ratio)
        f = exp(log(alp) + log(delta) - log(pi) + dg + bet*z + lk1 - 0.5_dp*log(d2pz2))
    end function pdf_nig_gen

    pure elemental function pdf_vg_sym(x, nu) result(f)
        ! Symmetric Variance Gamma, unit variance.  nu > 0.
        ! X = sqrt(G)*Z,  G ~ Gamma(1/nu, nu);  Var(X) = 1.
        ! f(x) = C * |x|^(1/nu-1/2) * K_{1/nu-1/2}(|x| * sqrt(2/nu))
        ! where C = sqrt(2/nu)^(1/nu+1/2) / (sqrt(pi)*Gamma(1/nu)*2^(1/nu-1/2)).
        real(dp), intent(in) :: x    ! argument
        real(dp), intent(in) :: nu   ! shape (> 0)
        real(dp) :: f, ord, s, ax, lc, lknu
        ord = 1.0_dp/nu - 0.5_dp
        s   = sqrt(2.0_dp / nu)
        ax  = abs(x)
        if (ax < 1.0e-300_dp) then
            f = 0.0_dp
            return
        end if
        lknu = bessel_k_nu(ord, ax * s)
        lc = (1.0_dp/nu + 0.5_dp)*log(s) &
             - 0.5_dp*log(pi) - log_gamma(1.0_dp/nu) &
             - (1.0_dp/nu - 0.5_dp)*log(2.0_dp)
        f  = exp(lc + (1.0_dp/nu - 0.5_dp)*log(ax) + lknu)
    end function pdf_vg_sym

    pure elemental function pdf_vg_gen(x, nu, rho) result(f)
        ! General (skewed) Variance Gamma, unit variance.
        ! nu > 0, rho = theta*sqrt(nu) in (-1,1).
        ! theta = rho/sqrt(nu),  sigma^2 = 1 - rho^2.
        ! X = theta*(G-1) + sigma*sqrt(G)*Z,  G ~ Gamma(1/nu, nu).
        ! Zero mean by construction; Var(X) = sigma^2 + theta^2*nu = 1.
        ! rho=0 recovers pdf_vg_sym(x, nu) exactly.
        real(dp), intent(in) :: x     ! argument
        real(dp), intent(in) :: nu    ! shape (> 0)
        real(dp), intent(in) :: rho   ! skewness in (-1,1)
        real(dp) :: f, ord, theta, sig2, z, r, lc, lknu
        ord   = 1.0_dp/nu - 0.5_dp
        theta = rho / sqrt(nu)
        sig2  = 1.0_dp - rho**2
        z     = x + theta             ! x - mu0, mu0 = -theta
        r     = sqrt(theta**2 + sig2 * z**2) / sig2
        if (r < 1.0e-300_dp) then
            f = 0.0_dp
            return
        end if
        lknu = bessel_k_nu(ord, r)
        lc = (1.0_dp/nu + 0.5_dp)*log(sqrt(2.0_dp/nu)) &
             - 0.5_dp*log(pi) - log_gamma(1.0_dp/nu) &
             - (1.0_dp/nu - 0.5_dp)*log(2.0_dp) &
             - log(sig2) &
             + (1.0_dp/nu - 0.5_dp)*log(0.5_dp * sig2)
        f = exp(lc + theta*z/sig2 + (1.0_dp/nu - 0.5_dp)*log(r) + lknu)
    end function pdf_vg_gen

    pure elemental function logpdf_std(x, dist_id, shape_val) result(logf)
        ! Log-density of a standardized innovation.
        real(dp), intent(in) :: x, shape_val
        integer,  intent(in) :: dist_id
        real(dp) :: logf
        logf = -nll_std_val(x, dist_id, shape_val)
    end function logpdf_std

    pure elemental function score_z_std(x, dist_id, shape_val) result(score)
        ! Score with respect to standardized residual z: d log f(z) / dz.
        real(dp), intent(in) :: x, shape_val
        integer,  intent(in) :: dist_id
        real(dp) :: score, nu, alp, lam, c, ax, a2px2, r, lk1, ratio, k1prime_over_k1
        select case (dist_id)
        case (dist_normal)
            score = -x
        case (dist_t)
            nu = shape_val
            score = -(nu + 1.0_dp) * x / (nu - 2.0_dp + x**2)
        case (dist_ged)
            nu = shape_val
            lam = ged_lambda(nu)
            ax = abs(x)
            if (ax == 0.0_dp) then
                score = 0.0_dp
            else
                score = -0.5_dp * nu * ax**(nu - 1.0_dp) * sign(1.0_dp, x) / lam**nu
            end if
        case (dist_logistic)
            c = pi / sqrt3
            score = -c * tanh(0.5_dp * c * x)
        case (dist_laplace)
            if (x == 0.0_dp) then
                score = 0.0_dp
            else
                score = -sqrt2 * sign(1.0_dp, x)
            end if
        case (dist_sech)
            score = -0.5_dp * pi * tanh(0.5_dp * pi * x)
        case (dist_nig_sym)
            alp = shape_val
            a2px2 = alp**2 + x**2
            r = alp * sqrt(a2px2)
            call bessel_k01(r, lk1, ratio)
            k1prime_over_k1 = -ratio - 1.0_dp / r
            score = x * (alp * k1prime_over_k1 / sqrt(a2px2) - 1.0_dp / a2px2)
        case default
            score = -x
        end select
    end function score_z_std

    pure elemental function score_shape_std(x, dist_id, shape_val) result(score)
        ! Score with respect to the constrained shape parameter: d log f(z;shape)/d shape.
        real(dp), intent(in) :: x, shape_val
        integer,  intent(in) :: dist_id
        real(dp) :: score, nu, u, h_step, shape_f, shape_b
        select case (dist_id)
        case (dist_t)
            nu = shape_val
            u = x**2 / (nu - 2.0_dp)
            score = 0.5_dp*digamma(0.5_dp*(nu + 1.0_dp)) &
                    - 0.5_dp*digamma(0.5_dp*nu) &
                    - 0.5_dp/(nu - 2.0_dp) &
                    - 0.5_dp*log(1.0_dp + u) &
                    + 0.5_dp*(nu + 1.0_dp)*u / ((nu - 2.0_dp)*(1.0_dp + u))
        case (dist_ged, dist_nig_sym)
            h_step = max(1.0e-5_dp * abs(shape_val), 1.0e-6_dp)
            shape_f = shape_val + h_step
            shape_b = max(shape_val - h_step, 1.0e-8_dp)
            score = (logpdf_std(x, dist_id, shape_f) - logpdf_std(x, dist_id, shape_b)) / &
                    (shape_f - shape_b)
        case default
            score = 0.0_dp
        end select
    end function score_shape_std

    ! ── Fitting: shape only (standardised residuals) ──────────────────────────

    subroutine fit_dist_std(x, n, dist_id, shape, loglik, converged, niter_out, max_iter_in, stats_start_in, &
                            start_exkurt_in, robust_start_in, start_tail_ratio_in, start_shape_in, xi_out, &
                            fixed_shape_in)
        ! Fit the shape parameter to standardised iid data x(1:n).
        ! mu=0 and sigma=1 are fixed; only shape is estimated (1-param BFGS).
        ! For no-shape distributions: NLL computed directly, converged=.true., shape=0.
        integer,  intent(in)  :: n, dist_id
        real(dp), intent(in)  :: x(n)
        real(dp), intent(out) :: shape, loglik
        logical,  intent(out) :: converged
        integer, optional, intent(out) :: niter_out
        integer, optional, intent(in) :: max_iter_in
        logical, optional, intent(in) :: stats_start_in
        real(dp), optional, intent(in) :: start_exkurt_in
        logical, optional, intent(in) :: robust_start_in
        real(dp), optional, intent(in) :: start_tail_ratio_in
        real(dp), optional, intent(in) :: start_shape_in
        real(dp), optional, intent(out) :: xi_out
        real(dp), optional, intent(in) :: fixed_shape_in
        integer,  parameter :: default_max_iter = 500
        real(dp), parameter :: gtol     = 1.0e-7_dp
        real(dp) :: p(2), fopt, exkurt, tail_ratio
        integer  :: niter, max_iter, t, np
        logical  :: stats_start, robust_start

        max_iter = default_max_iter
        if (present(max_iter_in)) max_iter = max_iter_in
        if (max_iter < 1) error stop "fit_dist_std: max_iter must be positive"
        stats_start = .false.
        if (present(stats_start_in)) stats_start = stats_start_in
        robust_start = .false.
        if (present(robust_start_in)) robust_start = robust_start_in

        if (allocated(fd_y)) deallocate(fd_y)
        allocate(fd_y(n))
        fd_y    = x
        fd_n    = n
        fd_dist = dist_id
        fd_mode = 1   ! shape-only mode
        fd_fixed_shape = .false.
        fd_fixed_shape_val = 0.0_dp

        if (present(xi_out)) xi_out = 0.0_dp
        if (present(fixed_shape_in)) then
            shape = fixed_shape_in
            converged = .true.
            if (present(niter_out)) niter_out = 0
            loglik = 0.0_dp
            do t = 1, n
                loglik = loglik - nll_std_val(x(t), dist_id, shape)
            end do
            return
        end if
        if (dist_npar_std(dist_id) == 0) then
            ! No shape parameter: compute NLL directly.
            shape     = 0.0_dp
            converged = .true.
            if (present(niter_out)) niter_out = 0
            loglik    = 0.0_dp
            do t = 1, n
                loglik = loglik - nll_std_val(x(t), dist_id, 0.0_dp)
            end do
            return
        end if

        ! 1-parameter BFGS on unconstrained shape.
        exkurt = 0.0_dp
        tail_ratio = 2.438_dp
        if (stats_start) then
            if (present(start_exkurt_in)) then
                exkurt = start_exkurt_in
            else
                exkurt = sample_excess_kurtosis(x)
            end if
        end if
        if (robust_start) then
            if (present(start_tail_ratio_in)) then
                tail_ratio = start_tail_ratio_in
            else
                tail_ratio = 2.438_dp
            end if
        end if
        select case (dist_id)
        case (dist_t)
            if (present(start_shape_in)) then
                p(1) = shape_inverse_transform(dist_id, start_shape_in)
            else if (robust_start) then
                p(1) = shape_inverse_transform(dist_id, shape_start_from_tail_ratio(dist_id, tail_ratio))
            else if (stats_start) then
                p(1) = shape_inverse_transform(dist_id, shape_start_from_exkurt(dist_id, exkurt))
            else
                p(1) = log((8.0_dp - 2.0_dp) / (100.0_dp - 8.0_dp))  ! nu = 8
            end if
        case (dist_ged)
            if (present(start_shape_in)) then
                p(1) = shape_inverse_transform(dist_id, start_shape_in)
            else if (robust_start) then
                p(1) = shape_inverse_transform(dist_id, shape_start_from_tail_ratio(dist_id, tail_ratio))
            else if (stats_start) then
                p(1) = shape_inverse_transform(dist_id, shape_start_from_exkurt(dist_id, exkurt))
            else
                p(1) = log(1.5_dp)                                   ! nu = 1.5
            end if
        case (dist_nig_sym)
            if (present(start_shape_in)) then
                p(1) = shape_inverse_transform(dist_id, start_shape_in)
            else if (robust_start) then
                if (present(start_exkurt_in)) then
                    p(1) = shape_inverse_transform(dist_id, shape_start_from_exkurt(dist_id, exkurt))
                else
                    p(1) = shape_inverse_transform(dist_id, shape_start_from_tail_ratio(dist_id, tail_ratio))
                end if
            else if (stats_start) then
                p(1) = shape_inverse_transform(dist_id, shape_start_from_exkurt(dist_id, exkurt))
            else
                p(1) = log((3.0_dp - 0.1_dp) / (20.0_dp - 3.0_dp))   ! alp = 3
            end if
        case (dist_nig_gen)
            if (present(start_shape_in)) then
                p(1) = shape_inverse_transform(dist_nig_sym, start_shape_in)
            else if (robust_start) then
                if (present(start_exkurt_in)) then
                    p(1) = shape_inverse_transform(dist_nig_sym, shape_start_from_exkurt(dist_nig_sym, exkurt))
                else
                    p(1) = shape_inverse_transform(dist_nig_sym, shape_start_from_tail_ratio(dist_nig_sym, tail_ratio))
                end if
            else if (stats_start) then
                p(1) = shape_inverse_transform(dist_nig_sym, shape_start_from_exkurt(dist_nig_sym, exkurt))
            else
                p(1) = log((3.0_dp - 0.1_dp) / (20.0_dp - 3.0_dp))   ! alp = 3
            end if
            p(2) = 0.0_dp                                             ! rho = beta/alp = 0
        case (dist_fs_skewt)
            if (present(start_shape_in)) then
                p(1) = shape_inverse_transform(dist_t, start_shape_in)
            else if (robust_start) then
                p(1) = shape_inverse_transform(dist_t, shape_start_from_tail_ratio(dist_t, tail_ratio))
            else if (stats_start) then
                p(1) = shape_inverse_transform(dist_t, shape_start_from_exkurt(dist_t, exkurt))
            else
                p(1) = log((8.0_dp - 2.0_dp) / (100.0_dp - 8.0_dp))  ! nu = 8
            end if
            p(2) = 0.0_dp                                             ! xi = 1
        case default
            p(1) = 0.0_dp
        end select

        np = max(1, dist_npar_std(dist_id))
        call bfgs_minimize(dist_obj, p, np, max_iter, gtol, fopt, niter, converged)
        if (present(niter_out)) niter_out = niter
        if (dist_id == dist_fs_skewt) then
            shape = shape_transform(dist_t, p(1))
            if (present(xi_out)) xi_out = exp(p(2))
        else if (dist_id == dist_nig_gen) then
            shape = shape_transform(dist_nig_sym, p(1))
            if (present(xi_out)) xi_out = rho_transform(p(2))
        else
            shape = shape_transform(dist_id, p(1))
        end if
        loglik = -fopt * real(n, dp)
    end subroutine fit_dist_std

    ! ── Fitting: mu + sigma + shape (raw iid data) ────────────────────────────

    subroutine fit_dist(y, n, dist_id, mu, sigma, shape, loglik, converged, niter_out, max_iter_in, stats_start_in, &
                        start_exkurt_in, robust_start_in, start_tail_ratio_in, start_shape_in, start_sigma_in, xi_out, &
                        fixed_shape_in)
        ! Fit distribution to raw iid data y(1:n).
        ! mu is fixed to the sample mean; sigma and shape estimated by BFGS.
        ! np = 1 (no-shape dists: just sigma) or 2 (t/ged/nig: sigma + shape).
        integer,  intent(in)  :: n, dist_id
        real(dp), intent(in)  :: y(n)
        real(dp), intent(out) :: mu, sigma, shape, loglik
        logical,  intent(out) :: converged
        integer, optional, intent(out) :: niter_out
        integer, optional, intent(in) :: max_iter_in
        logical, optional, intent(in) :: stats_start_in
        real(dp), optional, intent(in) :: start_exkurt_in
        logical, optional, intent(in) :: robust_start_in
        real(dp), optional, intent(in) :: start_tail_ratio_in
        real(dp), optional, intent(in) :: start_shape_in
        real(dp), optional, intent(in) :: start_sigma_in
        real(dp), optional, intent(out) :: xi_out
        real(dp), optional, intent(in) :: fixed_shape_in
        integer,  parameter :: default_max_iter = 500
        real(dp), parameter :: gtol     = 1.0e-7_dp
        real(dp) :: p(3), fopt, sample_sd, exkurt, tail_ratio
        integer  :: niter, max_iter, np
        logical  :: stats_start, robust_start

        max_iter = default_max_iter
        if (present(max_iter_in)) max_iter = max_iter_in
        if (max_iter < 1) error stop "fit_dist: max_iter must be positive"
        stats_start = .false.
        if (present(stats_start_in)) stats_start = stats_start_in
        robust_start = .false.
        if (present(robust_start_in)) robust_start = robust_start_in

        mu        = sum(y) / real(n, dp)
        sample_sd = sqrt(max(sum((y - mu)**2) / real(n, dp), 1.0e-20_dp))
        exkurt    = 0.0_dp
        tail_ratio = 2.438_dp
        if (stats_start) then
            if (present(start_exkurt_in)) then
                exkurt = start_exkurt_in
            else
                exkurt = sample_excess_kurtosis((y - mu) / sample_sd)
            end if
        end if
        if (robust_start) then
            if (present(start_tail_ratio_in)) then
                tail_ratio = start_tail_ratio_in
            else
                tail_ratio = 2.438_dp
            end if
        end if

        if (allocated(fd_y)) deallocate(fd_y)
        allocate(fd_y(n))
        fd_y    = y - mu
        fd_n    = n
        fd_dist = dist_id
        fd_mode = 0   ! sigma + shape mode
        fd_fixed_shape = .false.
        fd_fixed_shape_val = 0.0_dp
        if (present(fixed_shape_in)) then
            fd_fixed_shape = .true.
            fd_fixed_shape_val = fixed_shape_in
        end if
        if (present(xi_out)) xi_out = 0.0_dp

        if (present(start_sigma_in)) then
            p(1) = log(max(start_sigma_in*sample_sd, 1.0e-20_dp))
        else
            p(1) = log(sample_sd)
        end if
        if (present(fixed_shape_in)) then
            np = 1
            shape = fixed_shape_in
            call bfgs_minimize(dist_obj, p, np, max_iter, gtol, fopt, niter, converged)
            if (present(niter_out)) niter_out = niter
            sigma = exp(p(1))
            loglik = -fopt * real(n, dp)
            return
        end if
        select case (dist_id)
        case (dist_t)
            np   = 2
            if (present(start_shape_in)) then
                p(2) = shape_inverse_transform(dist_id, start_shape_in)
            else if (robust_start) then
                p(2) = shape_inverse_transform(dist_id, shape_start_from_tail_ratio(dist_id, tail_ratio))
            else if (stats_start) then
                p(2) = shape_inverse_transform(dist_id, shape_start_from_exkurt(dist_id, exkurt))
            else
                p(2) = log((8.0_dp - 2.0_dp) / (100.0_dp - 8.0_dp))  ! nu = 8
            end if
        case (dist_ged)
            np   = 2
            if (present(start_shape_in)) then
                p(2) = shape_inverse_transform(dist_id, start_shape_in)
            else if (robust_start) then
                p(2) = shape_inverse_transform(dist_id, shape_start_from_tail_ratio(dist_id, tail_ratio))
            else if (stats_start) then
                p(2) = shape_inverse_transform(dist_id, shape_start_from_exkurt(dist_id, exkurt))
            else
                p(2) = log(1.5_dp)                                      ! nu = 1.5
            end if
        case (dist_nig_sym)
            np   = 2
            if (present(start_shape_in)) then
                p(2) = shape_inverse_transform(dist_id, start_shape_in)
            else if (robust_start) then
                if (present(start_exkurt_in)) then
                    p(2) = shape_inverse_transform(dist_id, shape_start_from_exkurt(dist_id, exkurt))
                else
                    p(2) = shape_inverse_transform(dist_id, shape_start_from_tail_ratio(dist_id, tail_ratio))
                end if
            else if (stats_start) then
                p(2) = shape_inverse_transform(dist_id, shape_start_from_exkurt(dist_id, exkurt))
            else
                p(2) = log((3.0_dp - 0.1_dp) / (20.0_dp - 3.0_dp))   ! alp = 3
            end if
        case (dist_nig_gen)
            np = 3
            if (present(start_shape_in)) then
                p(2) = shape_inverse_transform(dist_nig_sym, start_shape_in)
            else if (robust_start) then
                if (present(start_exkurt_in)) then
                    p(2) = shape_inverse_transform(dist_nig_sym, shape_start_from_exkurt(dist_nig_sym, exkurt))
                else
                    p(2) = shape_inverse_transform(dist_nig_sym, shape_start_from_tail_ratio(dist_nig_sym, tail_ratio))
                end if
            else if (stats_start) then
                p(2) = shape_inverse_transform(dist_nig_sym, shape_start_from_exkurt(dist_nig_sym, exkurt))
            else
                p(2) = log((3.0_dp - 0.1_dp) / (20.0_dp - 3.0_dp))   ! alp = 3
            end if
            p(3) = 0.0_dp                                             ! rho = beta/alp = 0
        case (dist_fs_skewt)
            np = 3
            if (present(start_shape_in)) then
                p(2) = shape_inverse_transform(dist_t, start_shape_in)
            else if (robust_start) then
                p(2) = shape_inverse_transform(dist_t, shape_start_from_tail_ratio(dist_t, tail_ratio))
            else if (stats_start) then
                p(2) = shape_inverse_transform(dist_t, shape_start_from_exkurt(dist_t, exkurt))
            else
                p(2) = log((8.0_dp - 2.0_dp) / (100.0_dp - 8.0_dp))  ! nu = 8
            end if
            p(3) = 0.0_dp                                             ! xi = 1
        case default
            np   = 1
        end select

        call bfgs_minimize(dist_obj, p, np, max_iter, gtol, fopt, niter, converged)
        if (present(niter_out)) niter_out = niter
        sigma  = exp(p(1))
        shape  = 0.0_dp
        if (dist_id == dist_fs_skewt) then
            shape = shape_transform(dist_t, p(2))
            if (present(xi_out)) xi_out = exp(p(3))
        else if (dist_id == dist_nig_gen) then
            shape = shape_transform(dist_nig_sym, p(2))
            if (present(xi_out)) xi_out = rho_transform(p(3))
        else if (np == 2) then
            shape = shape_transform(dist_id, p(2))
        end if
        loglik = -fopt * real(n, dp)
    end subroutine fit_dist

    ! ── Private: BFGS callback ────────────────────────────────────────────────

    subroutine dist_obj(p, np, f, g)
        ! NLL/n with central-difference numerical gradient.
        ! Mode 0: p(1)=log(sigma), p(2)=unconstrained shape (if np=2)
        ! Mode 1: p(1)=unconstrained shape (sigma=1 fixed)
        integer,  intent(in)  :: np
        real(dp), intent(in)  :: p(np)
        real(dp), intent(out) :: f, g(np)
        real(dp), parameter   :: h_step = 1.0e-5_dp
        real(dp) :: pf(np), pb(np), ff, fb
        integer  :: j
        call dist_nll(p, np, f)
        do j = 1, np
            pf = p;  pf(j) = pf(j) + h_step
            pb = p;  pb(j) = pb(j) - h_step
            call dist_nll(pf, np, ff)
            call dist_nll(pb, np, fb)
            g(j) = (ff - fb) / (2.0_dp * h_step)
        end do
    end subroutine dist_obj

    subroutine dist_nll(p, np, f)
        ! Evaluate NLL/n from stored fd_y, fd_n, fd_dist, fd_mode.
        integer,  intent(in)  :: np
        real(dp), intent(in)  :: p(np)
        real(dp), intent(out) :: f
        real(dp) :: sigma, shape_val, xi_val, x, dens, s
        integer  :: t
        if (fd_mode == 0) then
            sigma     = exp(p(1))
            shape_val = 0.0_dp
            xi_val    = 1.0_dp
            if (fd_fixed_shape) then
                shape_val = fd_fixed_shape_val
            else if (fd_dist == dist_fs_skewt) then
                shape_val = shape_transform(dist_t, p(2))
                xi_val = exp(p(3))
            else if (fd_dist == dist_nig_gen) then
                shape_val = shape_transform(dist_nig_sym, p(2))
                xi_val = rho_transform(p(3))
            else if (np == 2) then
                shape_val = shape_transform(fd_dist, p(2))
            end if
            s = log(sigma)
        else
            sigma     = 1.0_dp
            xi_val    = 1.0_dp
            if (fd_fixed_shape) then
                shape_val = fd_fixed_shape_val
            else if (fd_dist == dist_fs_skewt) then
                shape_val = shape_transform(dist_t, p(1))
                xi_val = exp(p(2))
            else if (fd_dist == dist_nig_gen) then
                shape_val = shape_transform(dist_nig_sym, p(1))
                xi_val = rho_transform(p(2))
            else
                shape_val = shape_transform(fd_dist, p(1))
            end if
            s = 0.0_dp
        end if
        f = 0.0_dp
        do t = 1, fd_n
            x = fd_y(t) / sigma
            if (fd_dist == dist_fs_skewt) then
                dens = pdf_fs_skewt(x, shape_val, xi_val)
                f = f - log(max(dens, 1.0e-300_dp))
            else if (fd_dist == dist_nig_gen) then
                dens = pdf_nig_gen(x, shape_val, xi_val*shape_val)
                f = f - log(max(dens, 1.0e-300_dp))
            else
                f = f + nll_std_val(x, fd_dist, shape_val)
            end if
        end do
        f = s + f / real(fd_n, dp)
        if (f /= f .or. f > 1.0e29_dp) f = 1.0e30_dp   ! guard NaN/Inf
    end subroutine dist_nll

    ! ── Private: per-observation standardised NLL ─────────────────────────────

    pure elemental function nll_std_val(x, dist_id, shape_val) result(nll)
        ! NLL of a single standardised observation x under distribution dist_id.
        real(dp), intent(in) :: x, shape_val
        integer,  intent(in) :: dist_id
        real(dp) :: nll, nu, alp, lam, c, au, a2px2, lk1, ratio
        select case (dist_id)
        case (dist_normal)
            nll = log_sqrt_2pi + 0.5_dp*x**2
        case (dist_t)
            nu  = shape_val
            nll = log_gamma(0.5_dp*nu) - log_gamma(0.5_dp*(nu+1.0_dp)) &
                  + 0.5_dp*log(pi*(nu-2.0_dp)) &
                  + 0.5_dp*(nu+1.0_dp)*log(1.0_dp + x**2/(nu-2.0_dp))
        case (dist_ged)
            nu  = shape_val
            lam = ged_lambda(nu)
            nll = -log(nu) + log(2.0_dp) &
                  + 1.5_dp*log_gamma(1.0_dp/nu) - 0.5_dp*log_gamma(3.0_dp/nu) &
                  + 0.5_dp*(abs(x)/lam)**nu
        case (dist_logistic)
            c   = pi / sqrt3
            au  = c * abs(x)
            nll = -log(c) + au + 2.0_dp*log(1.0_dp + exp(-au))
        case (dist_laplace)
            nll = half_log2 + sqrt2*abs(x)
        case (dist_sech)
            au  = pi * abs(x) / 2.0_dp
            nll = au + log(1.0_dp + exp(-2.0_dp*au))   ! = log(2) + log(cosh(pi|x|/2))
        case (dist_nig_sym)
            alp   = shape_val
            a2px2 = alp**2 + x**2
            call bessel_k01(alp * sqrt(a2px2), lk1, ratio)
            nll = -2.0_dp*log(alp) + log(pi) - alp**2 + 0.5_dp*log(a2px2) - lk1
        case default
            nll = log_sqrt_2pi + 0.5_dp*x**2
        end select
    end function nll_std_val

    ! ── Private: shape parameter transforms ───────────────────────────────────

    pure elemental function shape_transform(dist_id, q) result(s)
        ! Unconstrained q -> constrained shape parameter.
        integer,  intent(in) :: dist_id
        real(dp), intent(in) :: q
        real(dp) :: s
        select case (dist_id)
        case (dist_t)
            s = 2.0_dp + 98.0_dp / (1.0_dp + exp(-q))        ! nu in (2, 100)
        case (dist_ged)
            s = exp(q)                                          ! nu > 0
        case (dist_nig_sym, dist_nig_gen)
            s = 0.1_dp + 19.9_dp / (1.0_dp + exp(-q))        ! alp in (0.1, 20)
        case default
            s = 0.0_dp
        end select
    end function shape_transform

    pure elemental function shape_inverse_transform(dist_id, s) result(q)
        ! Constrained shape parameter -> unconstrained optimiser parameter.
        integer, intent(in) :: dist_id
        real(dp), intent(in) :: s
        real(dp) :: q, sc

        select case (dist_id)
        case (dist_t)
            sc = min(max(s, 2.000001_dp), 99.999999_dp)
            q = log((sc - 2.0_dp) / (100.0_dp - sc))
        case (dist_ged)
            sc = max(s, 1.0e-8_dp)
            q = log(sc)
        case (dist_nig_sym, dist_nig_gen)
            sc = min(max(s, 0.100001_dp), 19.999999_dp)
            q = log((sc - 0.1_dp) / (20.0_dp - sc))
        case default
            q = 0.0_dp
        end select
    end function shape_inverse_transform

    pure elemental function rho_transform(q) result(rho)
        ! Unconstrained q -> NIG skewness ratio rho=beta/alpha in (-1,1).
        real(dp), intent(in) :: q
        real(dp) :: rho
        rho = 0.999_dp * tanh(q)
    end function rho_transform

    pure elemental function shape_start_from_exkurt(dist_id, exkurt) result(s)
        ! Conservative shape start based on sample excess kurtosis.
        integer, intent(in) :: dist_id
        real(dp), intent(in) :: exkurt
        real(dp) :: s, k

        k = max(exkurt, 0.0_dp)
        select case (dist_id)
        case (dist_t)
            if (k > 0.05_dp) then
                s = 4.0_dp + 6.0_dp / k
            else
                s = 50.0_dp
            end if
            s = min(max(s, 2.2_dp), 99.0_dp)
        case (dist_ged)
            if (exkurt >= 0.0_dp) then
                s = 2.0_dp / (1.0_dp + 0.35_dp*k)
            else
                s = 2.0_dp * (1.0_dp + min(abs(exkurt), 1.0_dp))
            end if
            s = min(max(s, 0.5_dp), 10.0_dp)
        case (dist_nig_sym, dist_nig_gen)
            if (k > 0.05_dp) then
                s = 3.0_dp / k
            else
                s = 10.0_dp
            end if
            s = min(max(s, 0.2_dp), 19.5_dp)
        case default
            s = 0.0_dp
        end select
    end function shape_start_from_exkurt

    pure elemental function shape_start_from_tail_ratio(dist_id, tail_ratio) result(s)
        ! Conservative shape start based on (q95-q05)/(q75-q25).
        integer, intent(in) :: dist_id
        real(dp), intent(in) :: tail_ratio
        real(dp) :: s, r, excess_ratio

        r = max(tail_ratio, 2.0_dp)
        excess_ratio = max(r - 2.438_dp, 0.0_dp)
        select case (dist_id)
        case (dist_t)
            if (excess_ratio > 0.05_dp) then
                s = 4.0_dp + 10.0_dp / excess_ratio
            else
                s = 50.0_dp
            end if
            s = min(max(s, 2.2_dp), 99.0_dp)
        case (dist_ged)
            if (r >= 2.438_dp) then
                s = 2.0_dp / (1.0_dp + 0.55_dp*excess_ratio)
            else
                s = 2.0_dp
            end if
            s = min(max(s, 0.5_dp), 10.0_dp)
        case (dist_nig_sym, dist_nig_gen)
            s = 3.0_dp
        case default
            s = 0.0_dp
        end select
    end function shape_start_from_tail_ratio

    pure elemental function dist_warm_shape_start(dist_id, best_simple_dist_id) result(s)
        ! KL-calibrated shape start from the best fitted no-shape distribution.
        integer, intent(in) :: dist_id, best_simple_dist_id
        real(dp) :: s

        select case (dist_id)
        case (dist_t)
            select case (best_simple_dist_id)
            case (dist_normal)
                s = 99.0_dp
            case (dist_logistic)
                s = 7.825601_dp
            case (dist_laplace)
                s = 3.803541_dp
            case (dist_sech)
                s = 5.300437_dp
            case default
                s = 8.0_dp
            end select
        case (dist_ged)
            select case (best_simple_dist_id)
            case (dist_normal)
                s = 2.0_dp
            case (dist_logistic)
                s = 1.463278_dp
            case (dist_laplace)
                s = 0.999993_dp
            case (dist_sech)
                s = 1.261982_dp
            case default
                s = 1.5_dp
            end select
        case (dist_nig_sym, dist_nig_gen)
            select case (best_simple_dist_id)
            case (dist_normal)
                s = 19.5_dp
            case (dist_logistic)
                s = 1.573270_dp
            case (dist_laplace)
                s = 0.795236_dp
            case (dist_sech)
                s = 1.149892_dp
            case default
                s = 3.0_dp
            end select
        case default
            s = 0.0_dp
        end select
    end function dist_warm_shape_start

    pure elemental function dist_warm_shape_start_scaled(dist_id, best_simple_dist_id) result(s)
        ! KL-calibrated shape start when target scale is also initialized.
        integer, intent(in) :: dist_id, best_simple_dist_id
        real(dp) :: s

        select case (dist_id)
        case (dist_t)
            select case (best_simple_dist_id)
            case (dist_normal)
                s = 99.0_dp
            case (dist_logistic)
                s = 7.763455_dp
            case (dist_laplace)
                s = 3.288721_dp
            case (dist_sech)
                s = 5.100085_dp
            case default
                s = 8.0_dp
            end select
        case (dist_ged)
            select case (best_simple_dist_id)
            case (dist_normal)
                s = 2.0_dp
            case (dist_logistic)
                s = 1.464218_dp
            case (dist_laplace)
                s = 0.999992_dp
            case (dist_sech)
                s = 1.264575_dp
            case default
                s = 1.5_dp
            end select
        case (dist_nig_sym, dist_nig_gen)
            select case (best_simple_dist_id)
            case (dist_normal)
                s = 19.5_dp
            case (dist_logistic)
                s = 1.573096_dp
            case (dist_laplace)
                s = 0.775065_dp
            case (dist_sech)
                s = 1.145753_dp
            case default
                s = 3.0_dp
            end select
        case default
            s = 0.0_dp
        end select
    end function dist_warm_shape_start_scaled

    pure elemental function dist_warm_sigma_start(dist_id, best_simple_dist_id) result(sigma)
        ! KL-calibrated scale start from the best fitted no-shape distribution.
        integer, intent(in) :: dist_id, best_simple_dist_id
        real(dp) :: sigma

        select case (dist_id)
        case (dist_t)
            select case (best_simple_dist_id)
            case (dist_normal)
                sigma = 1.000298_dp
            case (dist_logistic)
                sigma = 1.002222_dp
            case (dist_laplace)
                sigma = 1.086428_dp
            case (dist_sech)
                sigma = 1.013252_dp
            case default
                sigma = 1.0_dp
            end select
        case (dist_ged)
            select case (best_simple_dist_id)
            case (dist_normal)
                sigma = 1.000000_dp
            case (dist_logistic)
                sigma = 0.998553_dp
            case (dist_laplace)
                sigma = 1.000002_dp
            case (dist_sech)
                sigma = 0.996394_dp
            case default
                sigma = 1.0_dp
            end select
        case (dist_nig_sym, dist_nig_gen)
            select case (best_simple_dist_id)
            case (dist_normal)
                sigma = 1.000005_dp
            case (dist_logistic)
                sigma = 1.000057_dp
            case (dist_laplace)
                sigma = 1.016454_dp
            case (dist_sech)
                sigma = 1.002029_dp
            case default
                sigma = 1.0_dp
            end select
        case default
            sigma = 1.0_dp
        end select
    end function dist_warm_sigma_start

    pure function sample_excess_kurtosis(x) result(exkurt)
        ! Moment excess kurtosis using central moments with denominator n.
        real(dp), intent(in) :: x(:)
        real(dp) :: exkurt, mu, m2, m4
        integer :: n

        n = size(x)
        if (n < 1) then
            exkurt = 0.0_dp
            return
        end if
        mu = sum(x) / real(n, dp)
        m2 = sum((x - mu)**2) / real(n, dp)
        if (m2 <= 0.0_dp) then
            exkurt = 0.0_dp
        else
            m4 = sum((x - mu)**4) / real(n, dp)
            exkurt = m4 / m2**2 - 3.0_dp
        end if
    end function sample_excess_kurtosis

    integer function dist_id_from_name(dist_name) result(dist_id)
        ! Map a distribution name to its integer identifier.
        character(len=*), intent(in) :: dist_name

        select case (trim(uppercase_local(dist_name)))
        case ("NORMAL", "NORM", "GAUSSIAN")
            dist_id = dist_normal
        case ("T", "STUDENT", "STUDENT_T", "STUDENT-T")
            dist_id = dist_t
        case ("GED")
            dist_id = dist_ged
        case ("LOGISTIC")
            dist_id = dist_logistic
        case ("LAPLACE")
            dist_id = dist_laplace
        case ("SECH")
            dist_id = dist_sech
        case ("NIG_SYM", "SYM_NIG", "SYMMETRIC_NIG")
            dist_id = dist_nig_sym
        case ("NIG", "NIG_GEN", "GEN_NIG", "GENERAL_NIG")
            dist_id = dist_nig_gen
        case ("FS_SKEWT", "FS-SKEWT", "SKEWT")
            dist_id = dist_fs_skewt
        case default
            if (is_fixed_t_name(dist_name)) then
                dist_id = dist_t
            else
                dist_id = 0
            end if
        end select
    end function dist_id_from_name

    logical function dist_fixed_shape_from_name(dist_name, shape) result(has_fixed_shape)
        ! Parse fixed-shape distribution names such as T_6 or T_4P5.
        character(len=*), intent(in) :: dist_name
        real(dp), intent(out) :: shape
        character(len=len(dist_name)) :: name
        integer :: ios

        shape = 0.0_dp
        has_fixed_shape = .false.
        name = uppercase_local(adjustl(dist_name))
        if (.not. is_fixed_t_name(name)) return
        call read_real_token(name(3:), shape, ios)
        if (ios == 0 .and. shape > 2.0_dp) has_fixed_shape = .true.
    end function dist_fixed_shape_from_name

    logical function is_fixed_t_name(dist_name) result(is_fixed)
        ! Return true for names of the form T_# where # may use P for decimal point.
        character(len=*), intent(in) :: dist_name
        character(len=len(dist_name)) :: name
        real(dp) :: val
        integer :: ios

        name = uppercase_local(adjustl(dist_name))
        is_fixed = .false.
        if (len_trim(name) <= 2) return
        if (name(1:2) /= "T_") return
        call read_real_token(name(3:), val, ios)
        is_fixed = ios == 0 .and. val > 2.0_dp
    end function is_fixed_t_name

    subroutine read_real_token(token, value, ios)
        ! Read a real token, accepting P as a filename/header-safe decimal point.
        character(len=*), intent(in) :: token
        real(dp), intent(out) :: value
        integer, intent(out) :: ios
        character(len=len(token)) :: buf
        integer :: i

        buf = adjustl(token)
        do i = 1, len_trim(buf)
            if (buf(i:i) == "P") buf(i:i) = "."
        end do
        read(buf, *, iostat=ios) value
    end subroutine read_real_token

    pure function uppercase_local(s) result(out)
        ! Convert ASCII letters to uppercase without adding a module dependency.
        character(len=*), intent(in) :: s
        character(len=len(s)) :: out
        integer :: i, code

        out = s
        do i = 1, len(s)
            code = iachar(out(i:i))
            if (code >= iachar('a') .and. code <= iachar('z')) out(i:i) = achar(code - 32)
        end do
    end function uppercase_local

    ! ── Private: GED helpers ──────────────────────────────────────────────────

    pure elemental function ged_lambda(nu) result(lam)
        ! Scale for unit-variance GED: lam = 2^(-1/nu) * sqrt(Γ(1/nu)/Γ(3/nu))
        real(dp), intent(in) :: nu
        real(dp) :: lam
        lam = exp(-log(2.0_dp)/nu &
                  + 0.5_dp*(log_gamma(1.0_dp/nu) - log_gamma(3.0_dp/nu)))
    end function ged_lambda

end module distributions_mod

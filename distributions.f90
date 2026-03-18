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
!   dist_nig      (7)  Symmetric NIG, alp > 0    alp: tail/shape parameter
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
!
!   fit_dist(y, n, dist_id, mu, sigma, shape, loglik, converged)
!     Fit mu (= sample mean, fixed), sigma and shape to raw iid data y(1:n).
!     No-shape distributions: 1-parameter BFGS on sigma.
!     t, ged, nig: 2-parameter BFGS on (sigma, shape).
!
! Parameter counts for fit_dist_std (useful for AIC/BIC of shape comparisons):
!   dist_npar_std = [0, 1, 1, 0, 0, 0, 1]   (normal/logistic/laplace/sech = 0; t/ged/nig = 1)

module distributions_mod
    use kind_mod,       only: dp
    use math_const_mod, only: pi, two_pi, sqrt2, sqrt3, log2, half_log2, log_sqrt_2pi
    use bfgs_module,    only: bfgs_minimize
    use special_mod,    only: bessel_k01, bessel_k_nu
    implicit none
    private

    ! ── Distribution ID constants ─────────────────────────────────────────────

    integer, parameter, public :: dist_normal   = 1
    integer, parameter, public :: dist_t        = 2
    integer, parameter, public :: dist_ged      = 3
    integer, parameter, public :: dist_logistic = 4
    integer, parameter, public :: dist_laplace  = 5
    integer, parameter, public :: dist_sech     = 6
    integer, parameter, public :: dist_nig      = 7
    integer, parameter, public :: dist_count    = 7

    character(len=8), dimension(dist_count), parameter, public :: dist_names = &
        ["normal  ", "t       ", "ged     ", "logistic", &
         "laplace ", "sech    ", "nig     "]

    ! Number of shape parameters estimated by fit_dist_std (0 or 1 per distribution).
    integer, dimension(dist_count), parameter, public :: dist_npar_std = &
        [0, 1, 1, 0, 0, 0, 1]

    ! ── Module-level state for BFGS callbacks ─────────────────────────────────

    real(dp), allocatable, save :: fd_y(:)     ! stored (demeaned) data
    integer,               save :: fd_n    = 0
    integer,               save :: fd_dist = 0
    integer,               save :: fd_mode = 0 ! 0 = sigma+shape fit, 1 = shape-only fit

    ! ── Public interface ──────────────────────────────────────────────────────

    public :: pdf_normal, pdf_t, pdf_ged, pdf_logistic, pdf_laplace, pdf_sech, pdf_nig
    public :: pdf_nig_gen
    public :: pdf_vg_sym, pdf_vg_gen
    public :: fit_dist_std, fit_dist
    public :: ged_lambda

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

    ! ── Fitting: shape only (standardised residuals) ──────────────────────────

    subroutine fit_dist_std(x, n, dist_id, shape, loglik, converged)
        ! Fit the shape parameter to standardised iid data x(1:n).
        ! mu=0 and sigma=1 are fixed; only shape is estimated (1-param BFGS).
        ! For no-shape distributions: NLL computed directly, converged=.true., shape=0.
        integer,  intent(in)  :: n, dist_id
        real(dp), intent(in)  :: x(n)
        real(dp), intent(out) :: shape, loglik
        logical,  intent(out) :: converged
        integer,  parameter :: max_iter = 500
        real(dp), parameter :: gtol     = 1.0e-7_dp
        real(dp) :: p(1), fopt
        integer  :: niter, t

        if (allocated(fd_y)) deallocate(fd_y)
        allocate(fd_y(n))
        fd_y    = x
        fd_n    = n
        fd_dist = dist_id
        fd_mode = 1   ! shape-only mode

        if (dist_npar_std(dist_id) == 0) then
            ! No shape parameter: compute NLL directly.
            shape     = 0.0_dp
            converged = .true.
            loglik    = 0.0_dp
            do t = 1, n
                loglik = loglik - nll_std_val(x(t), dist_id, 0.0_dp)
            end do
            return
        end if

        ! 1-parameter BFGS on unconstrained shape.
        select case (dist_id)
        case (dist_t)
            p(1) = log((8.0_dp - 2.0_dp) / (100.0_dp - 8.0_dp))  ! nu = 8
        case (dist_ged)
            p(1) = log(1.5_dp)                                      ! nu = 1.5
        case (dist_nig)
            p(1) = log((3.0_dp - 0.1_dp) / (20.0_dp - 3.0_dp))   ! alp = 3
        case default
            p(1) = 0.0_dp
        end select

        call bfgs_minimize(dist_obj, p, 1, max_iter, gtol, fopt, niter, converged)
        shape  = shape_transform(dist_id, p(1))
        loglik = -fopt * real(n, dp)
    end subroutine fit_dist_std

    ! ── Fitting: mu + sigma + shape (raw iid data) ────────────────────────────

    subroutine fit_dist(y, n, dist_id, mu, sigma, shape, loglik, converged)
        ! Fit distribution to raw iid data y(1:n).
        ! mu is fixed to the sample mean; sigma and shape estimated by BFGS.
        ! np = 1 (no-shape dists: just sigma) or 2 (t/ged/nig: sigma + shape).
        integer,  intent(in)  :: n, dist_id
        real(dp), intent(in)  :: y(n)
        real(dp), intent(out) :: mu, sigma, shape, loglik
        logical,  intent(out) :: converged
        integer,  parameter :: max_iter = 500
        real(dp), parameter :: gtol     = 1.0e-7_dp
        real(dp) :: p(2), fopt, sample_sd
        integer  :: niter, np

        mu        = sum(y) / real(n, dp)
        sample_sd = sqrt(max(sum((y - mu)**2) / real(n, dp), 1.0e-20_dp))

        if (allocated(fd_y)) deallocate(fd_y)
        allocate(fd_y(n))
        fd_y    = y - mu
        fd_n    = n
        fd_dist = dist_id
        fd_mode = 0   ! sigma + shape mode

        p(1) = log(sample_sd)
        select case (dist_id)
        case (dist_t)
            np   = 2
            p(2) = log((8.0_dp - 2.0_dp) / (100.0_dp - 8.0_dp))  ! nu = 8
        case (dist_ged)
            np   = 2
            p(2) = log(1.5_dp)                                      ! nu = 1.5
        case (dist_nig)
            np   = 2
            p(2) = log((3.0_dp - 0.1_dp) / (20.0_dp - 3.0_dp))   ! alp = 3
        case default
            np   = 1
        end select

        call bfgs_minimize(dist_obj, p, np, max_iter, gtol, fopt, niter, converged)
        sigma  = exp(p(1))
        shape  = 0.0_dp
        if (np == 2) shape = shape_transform(dist_id, p(2))
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
        real(dp) :: sigma, shape_val, x, s
        integer  :: t
        if (fd_mode == 0) then
            sigma     = exp(p(1))
            shape_val = 0.0_dp
            if (np == 2) shape_val = shape_transform(fd_dist, p(2))
            s = log(sigma)
        else
            sigma     = 1.0_dp
            shape_val = shape_transform(fd_dist, p(1))
            s = 0.0_dp
        end if
        f = 0.0_dp
        do t = 1, fd_n
            x = fd_y(t) / sigma
            f = f + nll_std_val(x, fd_dist, shape_val)
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
        case (dist_nig)
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
        case (dist_nig)
            s = 0.1_dp + 19.9_dp / (1.0_dp + exp(-q))        ! alp in (0.1, 20)
        case default
            s = 0.0_dp
        end select
    end function shape_transform

    ! ── Private: GED helpers ──────────────────────────────────────────────────

    pure elemental function ged_lambda(nu) result(lam)
        ! Scale for unit-variance GED: lam = 2^(-1/nu) * sqrt(Γ(1/nu)/Γ(3/nu))
        real(dp), intent(in) :: nu
        real(dp) :: lam
        lam = exp(-log(2.0_dp)/nu &
                  + 0.5_dp*(log_gamma(1.0_dp/nu) - log_gamma(3.0_dp/nu)))
    end function ged_lambda

end module distributions_mod

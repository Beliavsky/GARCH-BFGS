module fgarch_mod
    ! Family GARCH: Hentschel (J. Fin. Econ. 39, 1995, pp. 71-104)
    !
    ! Variance recursion (normalised shock z_t = y_t / σ_t):
    !   h_t  = (σ_t^λ − 1) / λ                         [Box-Cox of σ]
    !   h_t  = ω + α · σ_{t-1}^ν · f(z_{t-1})^ν + β · h_{t-1}
    !   f(z) = |z − b| − c·(z − b)
    !
    ! Parameters p (7 for Normal, 8 for t):
    !   p(1) = ω  (free)
    !   p(2) = log α       → α = exp(p2) > 0
    !   p(3) = logit β     → β = sigmoid(p3) ∈ (0,1)
    !   p(4) = log λ       → λ = exp(p4) > 0
    !   p(5) = log ν       → ν = exp(p5) > 0
    !   p(6) = b           (free shift)
    !   p(7) = atanh c     → c = tanh(p7) ∈ (−1,1)
    !   p(8) = p_ν_t       → ν_t = 2 + 98/(1+exp(−p8))  [t only]
    !
    ! Gradient: analytical for p(1..7); 2-point FD for p(8) (t-df).

    use kind_mod, only: dp
    use math_const_mod, only: log_sqrt_2pi, sqrt2
    implicit none
    private

    public :: fg_dist_normal, fg_dist_t
    public :: fgarch_set_data, fgarch_set_dist, fgarch_np
    public :: fgarch_obj, fgarch_transform, fgarch_inv_transform
    public :: fgarch_vol_ann, fgarch_skew_kurt
    public :: fgarch22_obj, fgarch22_transform, fgarch22_inv_transform
    public :: fgarch22_vol_ann, fgarch22_skew_kurt
    public :: fgarch_twist_obj, fgarch_twist_transform, fgarch_twist_inv_transform
    public :: fgarch_twist_vol_ann, fgarch_twist_skew_kurt

    integer, parameter :: fg_dist_normal = 1
    integer, parameter :: fg_dist_t      = 2
    integer, parameter :: fg_np_base     = 7
    integer, parameter :: fg22_np        = 5
    integer, parameter :: fg_twist_np    = 5

    integer               :: fg_nobs = 0
    integer               :: fg_dist = fg_dist_normal
    real(dp), allocatable :: fg_y(:)

contains

    subroutine fgarch_set_data(y, n)
        real(dp), intent(in) :: y(n)
        integer,  intent(in) :: n
        fg_nobs = n
        if (allocated(fg_y)) deallocate(fg_y)
        allocate(fg_y(n))
        fg_y = y
    end subroutine

    subroutine fgarch_set_dist(idist)
        integer, intent(in) :: idist
        fg_dist = idist
    end subroutine

    integer function fgarch_np()
        if (fg_dist == fg_dist_t) then
            fgarch_np = fg_np_base + 1
        else
            fgarch_np = fg_np_base
        end if
    end function

    subroutine fgarch_transform(p, np, omega, alpha, beta, lam, nu, b, c, nu_t)
        integer,  intent(in)  :: np
        real(dp), intent(in)  :: p(np)
        real(dp), intent(out) :: omega, alpha, beta, lam, nu, b, c, nu_t
        omega = p(1)
        alpha = exp(p(2))
        beta  = 1.0_dp / (1.0_dp + exp(-p(3)))
        lam   = exp(p(4))
        nu    = exp(p(5))
        b     = p(6)
        c     = tanh(p(7))
        nu_t  = 0.0_dp
        if (np > fg_np_base) nu_t = 2.0_dp + 98.0_dp / (1.0_dp + exp(-p(8)))
    end subroutine

    subroutine fgarch_inv_transform(omega, alpha, beta, lam, nu, b, c, nu_t, p)
        real(dp), intent(in)  :: omega, alpha, beta, lam, nu, b, c, nu_t
        real(dp), intent(out) :: p(:)
        p(1) = omega
        p(2) = log(alpha)
        p(3) = log(beta / (1.0_dp - beta))
        p(4) = log(lam)
        p(5) = log(nu)
        p(6) = b
        p(7) = 0.5_dp * log((1.0_dp + c) / (1.0_dp - c))
        if (size(p) >= fg_np_base + 1) &
            p(8) = -log(98.0_dp / (nu_t - 2.0_dp) - 1.0_dp)
    end subroutine

    subroutine fgarch22_transform(p, omega, alpha, beta, b, c)
        real(dp), intent(in)  :: p(fg22_np)
        real(dp), intent(out) :: omega, alpha, beta, b, c
        omega = p(1)
        alpha = exp(p(2))
        beta  = 1.0_dp / (1.0_dp + exp(-p(3)))
        b     = p(4)
        c     = tanh(p(5))
    end subroutine

    subroutine fgarch22_inv_transform(omega, alpha, beta, b, c, p)
        real(dp), intent(in)  :: omega, alpha, beta, b, c
        real(dp), intent(out) :: p(fg22_np)
        p(1) = omega
        p(2) = log(alpha)
        p(3) = log(beta / (1.0_dp - beta))
        p(4) = b
        p(5) = 0.5_dp * log((1.0_dp + c) / (1.0_dp - c))
    end subroutine

    subroutine fgarch_twist_transform(p, omega, alpha, beta, theta, c)
        real(dp), intent(in)  :: p(fg_twist_np)
        real(dp), intent(out) :: omega, alpha, beta, theta, c
        real(dp) :: e2, e3, s, aa, moment, dmth, dmc
        omega = exp(p(1))
        theta = p(4)
        c = tanh(p(5))
        call fgarch_twist_moment(theta, c, moment, dmth, dmc)
        e2 = exp(p(2));  e3 = exp(p(3))
        s  = 1.0_dp + e2 + e3
        aa = e2 / s
        beta = e3 / s
        alpha = aa / moment
    end subroutine

    subroutine fgarch_twist_inv_transform(omega, alpha, beta, theta, c, p)
        real(dp), intent(in)  :: omega, alpha, beta, theta, c
        real(dp), intent(out) :: p(fg_twist_np)
        real(dp) :: aa, gam, moment, dmth, dmc
        call fgarch_twist_moment(theta, c, moment, dmth, dmc)
        aa = alpha * moment
        gam = 1.0_dp - aa - beta
        p(1) = log(omega)
        p(2) = log(aa / gam)
        p(3) = log(beta / gam)
        p(4) = theta
        p(5) = 0.5_dp * log((1.0_dp + c) / (1.0_dp - c))
    end subroutine

    ! NLL per obs only (used for FD on p(8)).
    subroutine fgarch_nll(p, np, f)
        integer,  intent(in)  :: np
        real(dp), intent(in)  :: p(np)
        real(dp), intent(out) :: f
        real(dp) :: omega, alpha, beta, lam, nu, b, c, nu_t
        real(dp) :: sigma, h, z, fz, fnu, h_new, slam, sigma_new, nll, lscore, s0
        real(dp), parameter :: pi = 3.14159265358979323846_dp
        integer  :: t
        call fgarch_transform(p, np, omega, alpha, beta, lam, nu, b, c, nu_t)
        s0    = sqrt(sum(fg_y**2) / fg_nobs)
        sigma = s0
        h     = (sigma**lam - 1.0_dp) / lam
        nll   = 0.0_dp
        do t = 1, fg_nobs
            z = fg_y(t) / sigma
            if (fg_dist == fg_dist_t) then
                lscore = log(sigma) &
                       + log_gamma(0.5_dp*nu_t) - log_gamma(0.5_dp*(nu_t+1.0_dp)) &
                       + 0.5_dp*log(pi*(nu_t-2.0_dp)) &
                       + 0.5_dp*(nu_t+1.0_dp)*log(1.0_dp + z**2/(nu_t-2.0_dp))
            else
                lscore = log_sqrt_2pi + log(sigma) + 0.5_dp*z**2
            end if
            nll  = nll + lscore
            fz   = abs(z - b) - c*(z - b)
            if (fz < 0.0_dp) fz = 0.0_dp
            fnu  = fz**nu
            h_new = omega + alpha * sigma**nu * fnu + beta * h
            slam  = 1.0_dp + lam * h_new
            if (slam <= 0.0_dp) then
                f = 1.0e30_dp
                return
            end if
            sigma_new = slam**(1.0_dp / lam)
            h     = h_new
            sigma = sigma_new
        end do
        f = nll / fg_nobs
    end subroutine

    ! NLL + analytical gradient for p(1..7).
    ! Tracks s = dh/dp and q = dσ/dp through the filter.
    subroutine fgarch_nll_grad(p, np, f, g7)
        integer,  intent(in)  :: np
        real(dp), intent(in)  :: p(np)
        real(dp), intent(out) :: f, g7(fg_np_base)
        real(dp) :: omega, alpha, beta, lam, nu, b, c, nu_t
        real(dp) :: sigma, h, z, fz, fnu, kappa, h_new, slam, sigma_new, nll, lscore, s0
        real(dp) :: Dfz, sign_zb, A_sigma, dsh, lam_dir, xi
        real(dp) :: direct(fg_np_base)
        real(dp) :: sv(fg_np_base), qv(fg_np_base), sv_n(fg_np_base), qv_n(fg_np_base)
        real(dp), parameter :: pi  = 3.14159265358979323846_dp
        real(dp), parameter :: eps = 1.0e-12_dp
        integer  :: t, k

        call fgarch_transform(p, np, omega, alpha, beta, lam, nu, b, c, nu_t)

        s0    = sqrt(sum(fg_y**2) / fg_nobs)
        sigma = s0
        h     = (sigma**lam - 1.0_dp) / lam

        ! sigma_0 = s0 is fixed → qv_0 = 0.
        ! h_0 = (s0^lam - 1)/lam depends on lam: dh_0/d(log_lam) = s0^lam*log(s0) - h_0.
        qv = 0.0_dp
        sv = 0.0_dp
        sv(4) = sigma**lam * log(sigma) - h    ! init sensitivity for log_lam

        nll = 0.0_dp
        g7  = 0.0_dp

        do t = 1, fg_nobs
            z = fg_y(t) / sigma

            ! NLL contribution and dl/dsigma (= xi).
            if (fg_dist == fg_dist_t) then
                lscore = log(sigma) &
                       + log_gamma(0.5_dp*nu_t) - log_gamma(0.5_dp*(nu_t+1.0_dp)) &
                       + 0.5_dp*log(pi*(nu_t-2.0_dp)) &
                       + 0.5_dp*(nu_t+1.0_dp)*log(1.0_dp + z**2/(nu_t-2.0_dp))
                xi = 1.0_dp/sigma &
                     - (nu_t+1.0_dp)*z**2 / (sigma*(nu_t-2.0_dp + z**2))
            else
                lscore = log_sqrt_2pi + log(sigma) + 0.5_dp*z**2
                xi = (1.0_dp - z**2) / sigma
            end if
            nll = nll + lscore

            ! Gradient from this observation: dl_t/dp = xi * qv.
            do k = 1, fg_np_base
                g7(k) = g7(k) + xi * qv(k)
            end do

            ! News impact at time t.
            if (z >= b) then
                sign_zb = 1.0_dp
            else
                sign_zb = -1.0_dp
            end if
            fz = abs(z - b) - c*(z - b)
            if (fz < 0.0_dp) fz = 0.0_dp

            if (fz > eps) then
                Dfz = nu * fz**(nu - 1.0_dp) * (sign_zb - c)
            else
                Dfz = 0.0_dp
            end if
            fnu   = fz**nu
            kappa = alpha * sigma**nu * fnu

            ! dh_{t+1}/dσ_t (holding h_t fixed).
            A_sigma = alpha * sigma**(nu - 1.0_dp) * (nu*fnu - Dfz*z)

            ! Direct contributions to dh_{t+1}/dp (holding σ_t, h_t fixed).
            direct(1) = 1.0_dp                                ! ω
            direct(2) = kappa                                 ! log α
            direct(3) = h * beta * (1.0_dp - beta)            ! logit β
            direct(4) = 0.0_dp                                ! log λ (none)
            if (fz > eps) then
                direct(5) = nu * kappa * (log(sigma) + log(fz))  ! log ν
                direct(6) = -alpha * sigma**nu * Dfz              ! b
                direct(7) = -alpha * sigma**nu &                  ! atanh c
                            * nu * fz**(nu-1.0_dp) * (z-b) * (1.0_dp - c**2)
            else
                direct(5) = 0.0_dp
                direct(6) = 0.0_dp
                direct(7) = 0.0_dp
            end if

            ! Update h and σ.
            h_new = omega + kappa + beta * h
            slam  = 1.0_dp + lam * h_new
            if (slam <= 0.0_dp) then
                f  = 1.0e30_dp
                g7 = 0.0_dp
                return
            end if
            sigma_new = slam**(1.0_dp / lam)
            dsh     = sigma_new**(1.0_dp - lam)          ! dsigma/dh = sigma^(1-lam)
            lam_dir = sigma_new * (lam * h_new / sigma_new**lam - log(sigma_new))

            ! Propagate sensitivities.
            do k = 1, fg_np_base
                sv_n(k) = direct(k) + beta * sv(k) + A_sigma * qv(k)
                qv_n(k) = dsh * sv_n(k)
            end do
            qv_n(4) = qv_n(4) + lam_dir    ! extra Box-Cox contribution to λ

            sv    = sv_n
            qv    = qv_n
            h     = h_new
            sigma = sigma_new
        end do

        f  = nll  / fg_nobs
        g7 = g7   / fg_nobs
    end subroutine

    ! Objective for BFGS: analytical gradient for p(1..7), 2-pt FD for p(8).
    subroutine fgarch_obj(p, np, f, g)
        integer,  intent(in)  :: np
        real(dp), intent(in)  :: p(np)
        real(dp), intent(out) :: f, g(np)
        real(dp) :: pp(np), fm, fp
        real(dp), parameter :: h_fd = 1.0e-5_dp

        call fgarch_nll_grad(p, np, f, g(1:fg_np_base))

        if (np > fg_np_base) then
            pp = p;  pp(fg_np_base+1) = p(fg_np_base+1) + h_fd
            call fgarch_nll(pp, np, fp)
            pp = p;  pp(fg_np_base+1) = p(fg_np_base+1) - h_fd
            call fgarch_nll(pp, np, fm)
            g(fg_np_base+1) = (fp - fm) / (2.0_dp * h_fd)
        end if
    end subroutine

    subroutine fgarch22_obj(p, np, f, g)
        integer,  intent(in)  :: np
        real(dp), intent(in)  :: p(np)
        real(dp), intent(out) :: f, g(np)
        real(dp) :: omega, alpha, beta, b, c
        real(dp) :: sigma, h, z, fz, fz2, kappa, h_new, slam, sigma_new, nll, lscore, s0
        real(dp) :: dfz2_dz, sign_zb, A_sigma, dsh, xi
        real(dp) :: direct(fg22_np), sv(fg22_np), qv(fg22_np), sv_n(fg22_np), qv_n(fg22_np)
        real(dp), parameter :: eps = 1.0e-12_dp
        integer :: t, k

        call fgarch22_transform(p(1:fg22_np), omega, alpha, beta, b, c)

        s0    = sqrt(sum(fg_y**2) / fg_nobs)
        sigma = s0
        h     = 0.5_dp * (sigma**2 - 1.0_dp)
        qv    = 0.0_dp
        sv    = 0.0_dp
        nll   = 0.0_dp
        g     = 0.0_dp

        do t = 1, fg_nobs
            z      = fg_y(t) / sigma
            lscore = log_sqrt_2pi + log(sigma) + 0.5_dp*z**2
            xi     = (1.0_dp - z**2) / sigma
            nll    = nll + lscore

            do k = 1, fg22_np
                g(k) = g(k) + xi * qv(k)
            end do

            if (z >= b) then
                sign_zb = 1.0_dp
            else
                sign_zb = -1.0_dp
            end if
            fz = abs(z - b) - c*(z - b)
            if (fz < 0.0_dp) fz = 0.0_dp

            if (fz > eps) then
                dfz2_dz = 2.0_dp * fz * (sign_zb - c)
            else
                dfz2_dz = 0.0_dp
            end if
            fz2   = fz**2
            kappa = alpha * sigma**2 * fz2

            A_sigma = alpha * sigma * (2.0_dp*fz2 - dfz2_dz*z)

            direct(1) = 1.0_dp
            direct(2) = kappa
            direct(3) = h * beta * (1.0_dp - beta)
            direct(4) = -alpha * sigma**2 * dfz2_dz
            if (fz > eps) then
                direct(5) = -alpha * sigma**2 * 2.0_dp * fz * (z - b) * (1.0_dp - c**2)
            else
                direct(5) = 0.0_dp
            end if

            h_new = omega + kappa + beta * h
            slam  = 1.0_dp + 2.0_dp * h_new
            if (slam <= 0.0_dp) then
                f = 1.0e30_dp
                g = 0.0_dp
                return
            end if
            sigma_new = sqrt(slam)
            dsh = 1.0_dp / sigma_new

            do k = 1, fg22_np
                sv_n(k) = direct(k) + beta * sv(k) + A_sigma * qv(k)
                qv_n(k) = dsh * sv_n(k)
            end do

            sv    = sv_n
            qv    = qv_n
            h     = h_new
            sigma = sigma_new
        end do

        f = nll / fg_nobs
        g = g / fg_nobs
    end subroutine

    ! Run filter; return annualised vol% = sqrt(252 * mean(σ_t^2)) * 100.
    subroutine fgarch_twist_obj(p, np, f, g)
        integer,  intent(in)  :: np
        real(dp), intent(in)  :: p(np)
        real(dp), intent(out) :: f, g(np)
        real(dp) :: omega, alpha, beta, theta, c
        real(dp) :: h, sqrth, z, x, q, q2, sign_x, kappa, nll, factor
        real(dp) :: moment, dmth, dmc, aa, D
        real(dp) :: dh_dom, dh_dal, dh_dbe, dh_dth, dh_dc
        real(dp) :: grad_om, grad_al, grad_be, grad_th, grad_c
        real(dp), parameter :: eps = 1.0e-12_dp
        integer :: t

        call fgarch_twist_transform(p(1:fg_twist_np), omega, alpha, beta, theta, c)
        call fgarch_twist_moment(theta, c, moment, dmth, dmc)

        D = 1.0_dp - alpha*moment - beta
        h = omega / D
        dh_dom = 1.0_dp / D
        dh_dal = h * moment / D
        dh_dbe = h / D
        dh_dth = h * alpha * dmth / D
        dh_dc  = h * alpha * dmc  / D

        nll = 0.0_dp
        grad_om = 0.0_dp;  grad_al = 0.0_dp;  grad_be = 0.0_dp
        grad_th = 0.0_dp;  grad_c  = 0.0_dp

        do t = 1, fg_nobs
            h = max(h, eps)
            sqrth = sqrt(h)
            z = fg_y(t) / sqrth
            nll = nll + log_sqrt_2pi + 0.5_dp * (log(h) + fg_y(t)**2 / h)
            factor = 1.0_dp/h - fg_y(t)**2 / h**2
            grad_om = grad_om + 0.5_dp * factor * dh_dom
            grad_al = grad_al + 0.5_dp * factor * dh_dal
            grad_be = grad_be + 0.5_dp * factor * dh_dbe
            grad_th = grad_th + 0.5_dp * factor * dh_dth
            grad_c  = grad_c  + 0.5_dp * factor * dh_dc

            x = z - theta
            sign_x = merge(1.0_dp, -1.0_dp, x >= 0.0_dp)
            q = abs(x) - c*x
            q2 = q*q
            kappa = beta + alpha * (q2 - q * (sign_x - c) * z)

            dh_dom = 1.0_dp                                  + kappa * dh_dom
            dh_dal = h * q2                                  + kappa * dh_dal
            dh_dbe = h                                       + kappa * dh_dbe
            dh_dth = -2.0_dp * alpha * h * q * (sign_x - c)  + kappa * dh_dth
            dh_dc  = -2.0_dp * alpha * h * q * x             + kappa * dh_dc
            h = omega + alpha * h * q2 + beta * h
        end do

        aa = alpha * moment
        g(1) = grad_om * omega
        g(2) = grad_al * alpha*(1.0_dp - aa) - grad_be * aa*beta
        g(3) = -grad_al * alpha*beta + grad_be * beta*(1.0_dp - beta)
        g(4) = grad_al * (-alpha*dmth/moment) + grad_th
        g(5) = (grad_al * (-alpha*dmc/moment) + grad_c) * (1.0_dp - c**2)

        f = nll / fg_nobs
        g = g / fg_nobs
    end subroutine

    subroutine fgarch_twist_moment(theta, c, moment, dmth, dmc)
        real(dp), intent(in)  :: theta, c
        real(dp), intent(out) :: moment, dmth, dmc
        real(dp) :: Phi, ph, A, B, epos, eneg, one_minus_c, one_plus_c

        Phi = 0.5_dp * (1.0_dp + erf(theta / sqrt2))
        ph  = exp(-0.5_dp * theta**2 - log_sqrt_2pi)

        A = (1.0_dp + theta**2) * (1.0_dp - Phi) - theta * ph
        B = (1.0_dp + theta**2) * Phi + theta * ph
        epos = ph - theta * (1.0_dp - Phi)
        eneg = -ph - theta * Phi

        one_minus_c = 1.0_dp - c
        one_plus_c  = 1.0_dp + c
        moment = one_minus_c**2 * A + one_plus_c**2 * B
        dmth   = -2.0_dp * (one_minus_c**2 * epos + one_plus_c**2 * eneg)
        dmc    = -2.0_dp * one_minus_c * A + 2.0_dp * one_plus_c * B
    end subroutine

    subroutine fgarch_vol_ann(p, np, vol_ann_pct)
        integer,  intent(in)  :: np
        real(dp), intent(in)  :: p(np)
        real(dp), intent(out) :: vol_ann_pct
        real(dp), parameter   :: trading_days = 252.0_dp
        real(dp) :: omega, alpha, beta, lam, nu, b, c, nu_t
        real(dp) :: sigma, h, z, fz, fnu, h_new, slam, sigma_new, s0, sum_var
        integer  :: t
        call fgarch_transform(p, np, omega, alpha, beta, lam, nu, b, c, nu_t)
        s0      = sqrt(sum(fg_y**2) / fg_nobs)
        sigma   = s0
        h       = (sigma**lam - 1.0_dp) / lam
        sum_var = 0.0_dp
        do t = 1, fg_nobs
            sum_var = sum_var + sigma**2
            z     = fg_y(t) / sigma
            fz    = abs(z - b) - c*(z - b)
            if (fz < 0.0_dp) fz = 0.0_dp
            fnu   = fz**nu
            h_new = omega + alpha * sigma**nu * fnu + beta * h
            slam  = 1.0_dp + lam * h_new
            if (slam <= 0.0_dp) slam = 1.0e-30_dp
            sigma_new = slam**(1.0_dp / lam)
            h     = h_new
            sigma = sigma_new
        end do
        vol_ann_pct = sqrt(trading_days * sum_var / fg_nobs) * 100.0_dp
    end subroutine

    subroutine fgarch22_vol_ann(p, np, vol_ann_pct)
        integer,  intent(in)  :: np
        real(dp), intent(in)  :: p(np)
        real(dp), intent(out) :: vol_ann_pct
        real(dp), parameter :: trading_days = 252.0_dp
        real(dp) :: omega, alpha, beta, b, c
        real(dp) :: sigma, h, z, fz, h_new, slam, sigma_new, s0, sum_var
        integer :: t
        call fgarch22_transform(p(1:fg22_np), omega, alpha, beta, b, c)
        s0      = sqrt(sum(fg_y**2) / fg_nobs)
        sigma   = s0
        h       = 0.5_dp * (sigma**2 - 1.0_dp)
        sum_var = 0.0_dp
        do t = 1, fg_nobs
            sum_var = sum_var + sigma**2
            z       = fg_y(t) / sigma
            fz      = abs(z - b) - c*(z - b)
            if (fz < 0.0_dp) fz = 0.0_dp
            h_new = omega + alpha * sigma**2 * fz**2 + beta * h
            slam  = 1.0_dp + 2.0_dp * h_new
            if (slam <= 0.0_dp) slam = 1.0e-30_dp
            sigma_new = sqrt(slam)
            h     = h_new
            sigma = sigma_new
        end do
        vol_ann_pct = sqrt(trading_days * sum_var / fg_nobs) * 100.0_dp
    end subroutine

    subroutine fgarch_twist_vol_ann(p, np, vol_ann_pct)
        integer,  intent(in)  :: np
        real(dp), intent(in)  :: p(np)
        real(dp), intent(out) :: vol_ann_pct
        real(dp), parameter :: trading_days = 252.0_dp
        real(dp) :: omega, alpha, beta, theta, c
        real(dp) :: h, z, x, q, moment, dmth, dmc, sum_var
        integer :: t

        call fgarch_twist_transform(p(1:fg_twist_np), omega, alpha, beta, theta, c)
        call fgarch_twist_moment(theta, c, moment, dmth, dmc)
        h = omega / (1.0_dp - alpha*moment - beta)
        sum_var = 0.0_dp
        do t = 1, fg_nobs
            h = max(h, 1.0e-12_dp)
            sum_var = sum_var + h
            z = fg_y(t) / sqrt(h)
            x = z - theta
            q = abs(x) - c*x
            h = omega + alpha * h * q**2 + beta * h
        end do
        vol_ann_pct = sqrt(trading_days * sum_var / fg_nobs) * 100.0_dp
    end subroutine

    ! Run filter; return skew and excess kurtosis of z_t = y_t / σ_t.
    subroutine fgarch_skew_kurt(p, np, skew, kurt)
        integer,  intent(in)  :: np
        real(dp), intent(in)  :: p(np)
        real(dp), intent(out) :: skew, kurt
        real(dp) :: omega, alpha, beta, lam, nu, b, c, nu_t
        real(dp) :: sigma, h, z, fz, fnu, h_new, slam, sigma_new, s0
        real(dp) :: zm, zv, zs, zk, dz, rn
        real(dp), allocatable :: zz(:)
        integer  :: t
        call fgarch_transform(p, np, omega, alpha, beta, lam, nu, b, c, nu_t)
        s0    = sqrt(sum(fg_y**2) / fg_nobs)
        sigma = s0
        h     = (sigma**lam - 1.0_dp) / lam
        allocate(zz(fg_nobs))
        do t = 1, fg_nobs
            z     = fg_y(t) / sigma
            zz(t) = z
            fz    = abs(z - b) - c*(z - b)
            if (fz < 0.0_dp) fz = 0.0_dp
            fnu   = fz**nu
            h_new = omega + alpha * sigma**nu * fnu + beta * h
            slam  = 1.0_dp + lam * h_new
            if (slam <= 0.0_dp) slam = 1.0e-30_dp
            sigma_new = slam**(1.0_dp / lam)
            h     = h_new
            sigma = sigma_new
        end do
        rn = real(fg_nobs, dp)
        zm = sum(zz) / rn
        zv = 0.0_dp;  zs = 0.0_dp;  zk = 0.0_dp
        do t = 1, fg_nobs
            dz = zz(t) - zm
            zv = zv + dz**2
            zs = zs + dz**3
            zk = zk + dz**4
        end do
        zv   = zv / rn
        skew = (zs / rn) / zv**1.5_dp
        kurt = (zk / rn) / zv**2 - 3.0_dp
        deallocate(zz)
    end subroutine

    subroutine fgarch22_skew_kurt(p, np, skew, kurt)
        integer,  intent(in)  :: np
        real(dp), intent(in)  :: p(np)
        real(dp), intent(out) :: skew, kurt
        real(dp) :: omega, alpha, beta, b, c
        real(dp) :: sigma, h, z, fz, h_new, slam, sigma_new, s0
        real(dp) :: zm, zv, zs, zk, dz, rn
        real(dp), allocatable :: zz(:)
        integer :: t
        call fgarch22_transform(p(1:fg22_np), omega, alpha, beta, b, c)
        s0    = sqrt(sum(fg_y**2) / fg_nobs)
        sigma = s0
        h     = 0.5_dp * (sigma**2 - 1.0_dp)
        allocate(zz(fg_nobs))
        do t = 1, fg_nobs
            z     = fg_y(t) / sigma
            zz(t) = z
            fz    = abs(z - b) - c*(z - b)
            if (fz < 0.0_dp) fz = 0.0_dp
            h_new = omega + alpha * sigma**2 * fz**2 + beta * h
            slam  = 1.0_dp + 2.0_dp * h_new
            if (slam <= 0.0_dp) slam = 1.0e-30_dp
            sigma_new = sqrt(slam)
            h     = h_new
            sigma = sigma_new
        end do
        rn = real(fg_nobs, dp)
        zm = sum(zz) / rn
        zv = 0.0_dp;  zs = 0.0_dp;  zk = 0.0_dp
        do t = 1, fg_nobs
            dz = zz(t) - zm
            zv = zv + dz**2
            zs = zs + dz**3
            zk = zk + dz**4
        end do
        zv   = zv / rn
        skew = (zs / rn) / zv**1.5_dp
        kurt = (zk / rn) / zv**2 - 3.0_dp
        deallocate(zz)
    end subroutine

    subroutine fgarch_twist_skew_kurt(p, np, skew, kurt)
        integer,  intent(in)  :: np
        real(dp), intent(in)  :: p(np)
        real(dp), intent(out) :: skew, kurt
        real(dp) :: omega, alpha, beta, theta, c
        real(dp) :: h, z, x, q, moment, dmth, dmc
        real(dp) :: zm, zv, zs, zk, dz, rn
        real(dp), allocatable :: zz(:)
        integer :: t

        call fgarch_twist_transform(p(1:fg_twist_np), omega, alpha, beta, theta, c)
        call fgarch_twist_moment(theta, c, moment, dmth, dmc)
        h = omega / (1.0_dp - alpha*moment - beta)
        allocate(zz(fg_nobs))
        do t = 1, fg_nobs
            h = max(h, 1.0e-12_dp)
            z = fg_y(t) / sqrt(h)
            zz(t) = z
            x = z - theta
            q = abs(x) - c*x
            h = omega + alpha * h * q**2 + beta * h
        end do
        rn = real(fg_nobs, dp)
        zm = sum(zz) / rn
        zv = 0.0_dp;  zs = 0.0_dp;  zk = 0.0_dp
        do t = 1, fg_nobs
            dz = zz(t) - zm
            zv = zv + dz**2
            zs = zs + dz**3
            zk = zk + dz**4
        end do
        zv   = zv / rn
        skew = (zs / rn) / zv**1.5_dp
        kurt = (zk / rn) / zv**2 - 3.0_dp
        deallocate(zz)
    end subroutine

end module fgarch_mod

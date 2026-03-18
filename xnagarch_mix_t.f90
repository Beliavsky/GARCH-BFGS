! NAGARCH with five noise distributions: Normal, Student-t, skewed-t, 2-component
! normal mixture (NM2), and 2-component Student-t mixture with shared dof (TM2).
! Branched from xnagarch_mix.f90.
!
! TM2 parameterisation (4 free params after mean=0, variance=1 constraints):
!   lambda = sigmoid(p5)          mixing weight of component 1
!   delta  = p6                   mean-separation: mu1 = delta*(1-lambda),
!                                                  mu2 = -delta*lambda
!   s      = sigmoid(p7)          within-variance split:
!                                  V = 1 - delta^2*lambda*(1-lambda)
!                                  sig_k^2 = s_k*V*(nu-2)/(nu*lambda_k)
!                                  (extra (nu-2)/nu vs NM2 for t-variance correction)
!   nu     = 2 + 98/(1+exp(-p8))  shared degrees of freedom
! Starting values p5=p6=p7=0, p8=p_t0 give lambda=0.5, delta=0, nu~8 -> t-mix collapses
! to symmetric equal-weight t with unit variance.

! ── NM2 module ───────────────────────────────────────────────────────────────
module nm2_nagarch_mod
    use kind_mod,       only: dp
    use math_const_mod, only: log_sqrt_2pi
    use nagarch_module, only: nagarch_transform
    implicit none
    private
    public :: nm2_set_data, nm2_obj, nm2_skew_kurt, nm2_nparams, nm2_unpack_dist

    integer, parameter :: nm2_nparams = 7   ! 4 NAGARCH + 3 mixture

    real(dp), allocatable, save :: nm2_obs(:)
    integer,               save :: nm2_nobs = 0

contains

    subroutine nm2_set_data(y, n)
        integer,  intent(in) :: n
        real(dp), intent(in) :: y(n)
        if (allocated(nm2_obs)) deallocate(nm2_obs)
        allocate(nm2_obs(n))
        nm2_obs  = y
        nm2_nobs = n
    end subroutine nm2_set_data

    ! Unpack unconstrained p(5:7) -> distribution moments (mu1,mu2,sig1,sig2,lambda).
    subroutine nm2_unpack_dist(p5, p6, p7, lam, mu1, mu2, sig1, sig2)
        real(dp), intent(in)  :: p5, p6, p7
        real(dp), intent(out) :: lam, mu1, mu2, sig1, sig2
        real(dp) :: delta, s, V
        real(dp), parameter :: eps = 1.0e-8_dp
        lam   = 1.0_dp / (1.0_dp + exp(-p5))
        delta = p6
        s     = 1.0_dp / (1.0_dp + exp(-p7))
        mu1   = delta * (1.0_dp - lam)
        mu2   = -delta * lam
        V     = max(1.0_dp - delta**2 * lam * (1.0_dp - lam), eps)
        sig1  = sqrt(s * V / max(lam, eps))
        sig2  = sqrt((1.0_dp - s) * V / max(1.0_dp - lam, eps))
    end subroutine nm2_unpack_dist

    ! NLL and d(log f)/dz for a 2-component normal mixture at standardised residual z.
    ! phi_k = exp(-u_k^2/2)/sig_k  (Gaussian kernel without 1/sqrt(2pi) factor)
    ! nll   = log_sqrt_2pi - log(lam*phi1 + (1-lam)*phi2)
    ! dl_dz = -(lam*phi1*u1/sig1 + (1-lam)*phi2*u2/sig2) / (lam*phi1 + (1-lam)*phi2)
    subroutine nm2_nll_score_z(z, p5, p6, p7, nll_z, dl_dz)
        real(dp), intent(in)  :: z, p5, p6, p7
        real(dp), intent(out) :: nll_z, dl_dz
        real(dp) :: lam, mu1, mu2, sig1, sig2
        real(dp) :: u1, u2, phi1, phi2, w1, w2, f_mix
        real(dp), parameter :: eps = 1.0e-300_dp
        call nm2_unpack_dist(p5, p6, p7, lam, mu1, mu2, sig1, sig2)
        u1    = (z - mu1) / sig1
        u2    = (z - mu2) / sig2
        phi1  = exp(-0.5_dp * u1**2) / sig1
        phi2  = exp(-0.5_dp * u2**2) / sig2
        w1    = lam * phi1
        w2    = (1.0_dp - lam) * phi2
        f_mix = max(w1 + w2, eps)
        nll_z = log_sqrt_2pi - log(f_mix)
        dl_dz = -(w1 * u1/sig1 + w2 * u2/sig2) / f_mix
    end subroutine nm2_nll_score_z

    ! NAGARCH-NM2 objective: NLL/n and gradient w.r.t. unconstrained p(1:7).
    ! p(1:4): NAGARCH structural (via nagarch_transform / nagarch_inv_transform)
    ! p(5:7): NM2 distribution params; gradient via central differences.
    subroutine nm2_obj(p, np, f, g)
        integer,  intent(in)  :: np
        real(dp), intent(in)  :: p(np)
        real(dp), intent(out) :: f, g(np)

        real(dp) :: omega, alpha, beta, theta, s2, aa, D
        real(dp) :: h, sqrth, r, z, nll_z, dl_dz, factor, kappa
        real(dp) :: dh_dom, dh_dal, dh_dbe, dh_dth
        real(dp) :: gs(4)        ! structural gradient accumulators
        real(dp) :: gm(3)        ! mixture gradient accumulators
        real(dp) :: nll_p, nll_m, dl_dummy
        real(dp), parameter :: h_fd  = 1.0e-5_dp
        real(dp), parameter :: eps_h = 1.0e-12_dp
        integer  :: t

        call nagarch_transform(p(1:4), omega, alpha, beta, theta)
        s2 = 1.0_dp + theta**2
        D  = max(1.0_dp - alpha*s2 - beta, 1.0e-8_dp)
        h  = omega / D

        ! Stationarity-based initial derivatives
        dh_dom = 1.0_dp / D
        dh_dal = h * s2 / D
        dh_dbe = h / D
        dh_dth = 2.0_dp * theta * alpha * h / D

        f   = 0.0_dp
        gs  = 0.0_dp
        gm  = 0.0_dp

        do t = 1, nm2_nobs
            sqrth = sqrt(max(h, eps_h))
            z     = nm2_obs(t) / sqrth
            call nm2_nll_score_z(z, p(5), p(6), p(7), nll_z, dl_dz)
            factor = (1.0_dp + dl_dz * z) / (2.0_dp * max(h, eps_h))
            f      = f      + nll_z + 0.5_dp * log(max(h, eps_h))
            gs(1)  = gs(1)  + factor * dh_dom
            gs(2)  = gs(2)  + factor * dh_dal
            gs(3)  = gs(3)  + factor * dh_dbe
            gs(4)  = gs(4)  + factor * dh_dth

            ! Central-difference gradient for NM2 params
            call nm2_nll_score_z(z, p(5)+h_fd, p(6),      p(7),      nll_p, dl_dummy)
            call nm2_nll_score_z(z, p(5)-h_fd, p(6),      p(7),      nll_m, dl_dummy)
            gm(1) = gm(1) + (nll_p - nll_m) / (2.0_dp * h_fd)

            call nm2_nll_score_z(z, p(5),      p(6)+h_fd, p(7),      nll_p, dl_dummy)
            call nm2_nll_score_z(z, p(5),      p(6)-h_fd, p(7),      nll_m, dl_dummy)
            gm(2) = gm(2) + (nll_p - nll_m) / (2.0_dp * h_fd)

            call nm2_nll_score_z(z, p(5),      p(6),      p(7)+h_fd, nll_p, dl_dummy)
            call nm2_nll_score_z(z, p(5),      p(6),      p(7)-h_fd, nll_m, dl_dummy)
            gm(3) = gm(3) + (nll_p - nll_m) / (2.0_dp * h_fd)

            ! Advance NAGARCH recurrence
            r      = nm2_obs(t) - theta * sqrth
            kappa  = beta - alpha * theta * r / sqrth
            dh_dom = 1.0_dp               + kappa * dh_dom
            dh_dal = r**2                 + kappa * dh_dal
            dh_dbe = h                    + kappa * dh_dbe
            dh_dth = -2.0_dp*alpha*r*sqrth + kappa * dh_dth
            h      = omega + alpha*r**2 + beta*h
        end do

        f   = f   / nm2_nobs
        gs  = gs  / nm2_nobs
        gm  = gm  / nm2_nobs

        ! Chain rule: structural -> unconstrained (same as proc_nagarch in garch_flex_mod)
        aa   = alpha * s2
        g(1) =  gs(1) * omega
        g(2) =  gs(2) * alpha*(1.0_dp - aa) - gs(3) * aa*beta
        g(3) =  gs(2) * (-alpha*beta)        + gs(3) * beta*(1.0_dp - beta)
        g(4) =  gs(2) * (-2.0_dp*theta*alpha/s2) + gs(4)
        g(5) =  gm(1)
        g(6) =  gm(2)
        g(7) =  gm(3)
    end subroutine nm2_obj

    ! Standardised residuals skewness and excess kurtosis for NAGARCH-NM2.
    subroutine nm2_skew_kurt(p, skew, kurt)
        real(dp), intent(in)  :: p(nm2_nparams)
        real(dp), intent(out) :: skew, kurt
        real(dp) :: omega, alpha, beta, theta, h, sqrth, r
        real(dp), allocatable :: zz(:)
        real(dp) :: zm, zv, zs, zk, dz, dz2, rn
        integer  :: t
        call nagarch_transform(p(1:4), omega, alpha, beta, theta)
        h = omega / max(1.0_dp - alpha*(1.0_dp+theta**2) - beta, 1.0e-8_dp)
        allocate(zz(nm2_nobs))
        do t = 1, nm2_nobs
            sqrth  = sqrt(max(h, 1.0e-12_dp))
            zz(t)  = nm2_obs(t) / sqrth
            r      = nm2_obs(t) - theta * sqrth
            h      = omega + alpha*r**2 + beta*h
        end do
        rn = real(nm2_nobs, dp)
        zm = sum(zz) / rn
        zv = 0.0_dp;  zs = 0.0_dp;  zk = 0.0_dp
        do t = 1, nm2_nobs
            dz  = zz(t) - zm
            dz2 = dz**2
            zv  = zv + dz2
            zs  = zs + dz*dz2
            zk  = zk + dz2**2
        end do
        zv   = zv / rn
        skew = (zs/rn) / zv**1.5_dp
        kurt = (zk/rn) / zv**2 - 3.0_dp
        deallocate(zz)
    end subroutine nm2_skew_kurt

end module nm2_nagarch_mod


! ── TM2 module ───────────────────────────────────────────────────────────────
module tm2_nagarch_mod
    use kind_mod,       only: dp
    use nagarch_module, only: nagarch_transform
    implicit none
    private
    public :: tm2_set_data, tm2_obj, tm2_skew_kurt, tm2_nparams, tm2_unpack_dist

    integer, parameter :: tm2_nparams = 8   ! 4 NAGARCH + 4 mixture

    real(dp), allocatable, save :: tm2_obs(:)
    integer,               save :: tm2_nobs = 0

contains

    subroutine tm2_set_data(y, n)
        integer,  intent(in) :: n
        real(dp), intent(in) :: y(n)
        if (allocated(tm2_obs)) deallocate(tm2_obs)
        allocate(tm2_obs(n))
        tm2_obs  = y
        tm2_nobs = n
    end subroutine tm2_set_data

    ! Unpack unconstrained p(5:8) -> distribution params.
    ! sig_k^2 scaled by (nu-2)/nu so the t-mixture has unit variance.
    subroutine tm2_unpack_dist(p5, p6, p7, p8, lam, mu1, mu2, sig1, sig2, nu)
        real(dp), intent(in)  :: p5, p6, p7, p8
        real(dp), intent(out) :: lam, mu1, mu2, sig1, sig2, nu
        real(dp) :: delta, s, V, scale
        real(dp), parameter :: eps = 1.0e-8_dp
        lam   = 1.0_dp / (1.0_dp + exp(-p5))
        delta = p6
        s     = 1.0_dp / (1.0_dp + exp(-p7))
        nu    = 2.0_dp + 98.0_dp / (1.0_dp + exp(-p8))
        mu1   = delta * (1.0_dp - lam)
        mu2   = -delta * lam
        V     = max(1.0_dp - delta**2 * lam * (1.0_dp - lam), eps)
        scale = (nu - 2.0_dp) / nu
        sig1  = sqrt(s * V * scale / max(lam, eps))
        sig2  = sqrt((1.0_dp - s) * V * scale / max(1.0_dp - lam, eps))
    end subroutine tm2_unpack_dist

    ! NLL and d(log f)/dz for a 2-component t mixture at standardised residual z.
    ! Each component is t(mu_k, sig_k, nu); shared nu.
    ! g_k = (1 + u_k^2/nu)^(-(nu+1)/2) / sig_k  (t kernel, nu-dependent normalisation
    !        absorbed into log_ct computed once per call)
    ! nll_z = -(log_ct + log(lam*g1 + (1-lam)*g2))
    ! dl_dz = -(lam*g1*(nu+1)*u1/sig1/(nu+u1^2) + (1-lam)*g2*(nu+1)*u2/sig2/(nu+u2^2))
    !          / f_mix
    subroutine tm2_nll_score_z(z, p5, p6, p7, p8, nll_z, dl_dz)
        real(dp), intent(in)  :: z, p5, p6, p7, p8
        real(dp), intent(out) :: nll_z, dl_dz
        real(dp) :: lam, mu1, mu2, sig1, sig2, nu
        real(dp) :: u1, u2, g1, g2, w1, w2, f_mix, log_ct
        real(dp), parameter :: eps = 1.0e-300_dp
        call tm2_unpack_dist(p5, p6, p7, p8, lam, mu1, mu2, sig1, sig2, nu)
        log_ct = log_gamma(0.5_dp*(nu+1.0_dp)) - log_gamma(0.5_dp*nu) &
                 - 0.5_dp*log(acos(-1.0_dp)*nu)
        u1    = (z - mu1) / sig1
        u2    = (z - mu2) / sig2
        g1    = (1.0_dp + u1**2/nu)**(-(nu+1.0_dp)*0.5_dp) / sig1
        g2    = (1.0_dp + u2**2/nu)**(-(nu+1.0_dp)*0.5_dp) / sig2
        w1    = lam * g1
        w2    = (1.0_dp - lam) * g2
        f_mix = max(w1 + w2, eps)
        nll_z = -(log_ct + log(f_mix))
        dl_dz = -(w1*(nu+1.0_dp)*u1/sig1/(nu+u1**2) &
                + w2*(nu+1.0_dp)*u2/sig2/(nu+u2**2)) / f_mix
    end subroutine tm2_nll_score_z

    ! NAGARCH-TM2 objective: NLL/n and gradient w.r.t. unconstrained p(1:8).
    ! p(1:4): NAGARCH structural params.
    ! p(5:8): TM2 distribution params; gradient via central differences.
    subroutine tm2_obj(p, np, f, g)
        integer,  intent(in)  :: np
        real(dp), intent(in)  :: p(np)
        real(dp), intent(out) :: f, g(np)

        real(dp) :: omega, alpha, beta, theta, s2, aa, D
        real(dp) :: h, sqrth, r, z, nll_z, dl_dz, factor, kappa
        real(dp) :: dh_dom, dh_dal, dh_dbe, dh_dth
        real(dp) :: gs(4)        ! structural gradient accumulators
        real(dp) :: gm(4)        ! mixture gradient accumulators
        real(dp) :: nll_p, nll_m, dl_dummy
        real(dp), parameter :: h_fd  = 1.0e-5_dp
        real(dp), parameter :: eps_h = 1.0e-12_dp
        integer  :: t

        call nagarch_transform(p(1:4), omega, alpha, beta, theta)
        s2 = 1.0_dp + theta**2
        D  = max(1.0_dp - alpha*s2 - beta, 1.0e-8_dp)
        h  = omega / D

        ! Stationarity-based initial derivatives
        dh_dom = 1.0_dp / D
        dh_dal = h * s2 / D
        dh_dbe = h / D
        dh_dth = 2.0_dp * theta * alpha * h / D

        f   = 0.0_dp
        gs  = 0.0_dp
        gm  = 0.0_dp

        do t = 1, tm2_nobs
            sqrth = sqrt(max(h, eps_h))
            z     = tm2_obs(t) / sqrth
            call tm2_nll_score_z(z, p(5), p(6), p(7), p(8), nll_z, dl_dz)
            factor = (1.0_dp + dl_dz * z) / (2.0_dp * max(h, eps_h))
            f      = f      + nll_z + 0.5_dp * log(max(h, eps_h))
            gs(1)  = gs(1)  + factor * dh_dom
            gs(2)  = gs(2)  + factor * dh_dal
            gs(3)  = gs(3)  + factor * dh_dbe
            gs(4)  = gs(4)  + factor * dh_dth

            ! Central-difference gradient for TM2 params
            call tm2_nll_score_z(z, p(5)+h_fd, p(6),      p(7),      p(8),      nll_p, dl_dummy)
            call tm2_nll_score_z(z, p(5)-h_fd, p(6),      p(7),      p(8),      nll_m, dl_dummy)
            gm(1) = gm(1) + (nll_p - nll_m) / (2.0_dp * h_fd)

            call tm2_nll_score_z(z, p(5),      p(6)+h_fd, p(7),      p(8),      nll_p, dl_dummy)
            call tm2_nll_score_z(z, p(5),      p(6)-h_fd, p(7),      p(8),      nll_m, dl_dummy)
            gm(2) = gm(2) + (nll_p - nll_m) / (2.0_dp * h_fd)

            call tm2_nll_score_z(z, p(5),      p(6),      p(7)+h_fd, p(8),      nll_p, dl_dummy)
            call tm2_nll_score_z(z, p(5),      p(6),      p(7)-h_fd, p(8),      nll_m, dl_dummy)
            gm(3) = gm(3) + (nll_p - nll_m) / (2.0_dp * h_fd)

            call tm2_nll_score_z(z, p(5),      p(6),      p(7),      p(8)+h_fd, nll_p, dl_dummy)
            call tm2_nll_score_z(z, p(5),      p(6),      p(7),      p(8)-h_fd, nll_m, dl_dummy)
            gm(4) = gm(4) + (nll_p - nll_m) / (2.0_dp * h_fd)

            ! Advance NAGARCH recurrence
            r      = tm2_obs(t) - theta * sqrth
            kappa  = beta - alpha * theta * r / sqrth
            dh_dom = 1.0_dp                + kappa * dh_dom
            dh_dal = r**2                  + kappa * dh_dal
            dh_dbe = h                     + kappa * dh_dbe
            dh_dth = -2.0_dp*alpha*r*sqrth + kappa * dh_dth
            h      = omega + alpha*r**2 + beta*h
        end do

        f   = f   / tm2_nobs
        gs  = gs  / tm2_nobs
        gm  = gm  / tm2_nobs

        ! Chain rule: structural -> unconstrained (same as proc_nagarch in garch_flex_mod)
        aa   = alpha * s2
        g(1) =  gs(1) * omega
        g(2) =  gs(2) * alpha*(1.0_dp - aa) - gs(3) * aa*beta
        g(3) =  gs(2) * (-alpha*beta)        + gs(3) * beta*(1.0_dp - beta)
        g(4) =  gs(2) * (-2.0_dp*theta*alpha/s2) + gs(4)
        g(5) =  gm(1)
        g(6) =  gm(2)
        g(7) =  gm(3)
        g(8) =  gm(4)
    end subroutine tm2_obj

    ! Standardised residuals skewness and excess kurtosis for NAGARCH-TM2.
    subroutine tm2_skew_kurt(p, skew, kurt)
        real(dp), intent(in)  :: p(tm2_nparams)
        real(dp), intent(out) :: skew, kurt
        real(dp) :: omega, alpha, beta, theta, h, sqrth, r
        real(dp), allocatable :: zz(:)
        real(dp) :: zm, zv, zs, zk, dz, dz2, rn
        integer  :: t
        call nagarch_transform(p(1:4), omega, alpha, beta, theta)
        h = omega / max(1.0_dp - alpha*(1.0_dp+theta**2) - beta, 1.0e-8_dp)
        allocate(zz(tm2_nobs))
        do t = 1, tm2_nobs
            sqrth  = sqrt(max(h, 1.0e-12_dp))
            zz(t)  = tm2_obs(t) / sqrth
            r      = tm2_obs(t) - theta * sqrth
            h      = omega + alpha*r**2 + beta*h
        end do
        rn = real(tm2_nobs, dp)
        zm = sum(zz) / rn
        zv = 0.0_dp;  zs = 0.0_dp;  zk = 0.0_dp
        do t = 1, tm2_nobs
            dz  = zz(t) - zm
            dz2 = dz**2
            zv  = zv + dz2
            zs  = zs + dz*dz2
            zk  = zk + dz2**2
        end do
        zv   = zv / rn
        skew = (zs/rn) / zv**1.5_dp
        kurt = (zk/rn) / zv**2 - 3.0_dp
        deallocate(zz)
    end subroutine tm2_skew_kurt

end module tm2_nagarch_mod


! ── Main program ─────────────────────────────────────────────────────────────
program xnagarch_mix_t
! Fit NAGARCH with Normal, Student-t, skewed-t, 2-component normal mixture (NM2),
! and 2-component Student-t mixture with shared dof (TM2).

use kind_mod,         only: dp
use csv_mod,          only: read_price_csv
use garch_flex_mod,   only: flex_set_data, flex_set_types, flex_np, flex_obj, flex_skew_kurt, &
                              proc_nagarch, dist_normal, dist_t, dist_skew_t
use nagarch_module,   only: nagarch_transform, nagarch_inv_transform
use nm2_nagarch_mod,  only: nm2_set_data, nm2_obj, nm2_skew_kurt, nm2_nparams, nm2_unpack_dist
use tm2_nagarch_mod,  only: tm2_set_data, tm2_obj, tm2_skew_kurt, tm2_nparams, tm2_unpack_dist
use bfgs_module,      only: bfgs_minimize
use stats_mod,        only: mean, sd
use rank_mod,         only: rank_desc, rank_asc
implicit none

integer,  parameter :: n_dist       = 5          ! Normal, t, skew-t, NM2, TM2
real(dp), parameter :: trading_days = 252.0_dp
integer,  parameter :: max_iter     = 500
real(dp), parameter :: gtol         = 1.0e-6_dp
integer,  parameter :: nret         = 10**6
real(dp), parameter :: p_t0         = -2.729_dp  ! starting logit for nu ~ 8

character(len=9), parameter :: dist_names(n_dist) = &
    ["Normal   ", "Student-t", "Skewed-t ", "NM2      ", "TM2      "]

! data
integer,           allocatable :: dates(:)
character(len=32), allocatable :: col_names(:)
real(dp),          allocatable :: prices(:,:)
real(dp),          allocatable :: raw_ret(:), ret(:)
integer :: nobs, nprices, nall, ncols, icol, i1

real(dp) :: ret_mean, ret_std

! optimisation
real(dp), allocatable :: p(:), p0(:)
real(dp) :: fopt
integer  :: np, niter, id
logical  :: converged
real(dp) :: t_start, t_end

! per-model results
real(dp) :: omegas(n_dist), alphas(n_dist), betas(n_dist), thetas(n_dist)
real(dp) :: d1s(n_dist), d2s(n_dist), d3s(n_dist), d4s(n_dist)  ! distribution display params
real(dp) :: mu1s(n_dist), mu2s(n_dist)                           ! component means (0 for non-mixture)
real(dp) :: vol_anns(n_dist), logls(n_dist), aics(n_dist), bics(n_dist)
real(dp) :: skews(n_dist), kurts(n_dist)
integer  :: niters(n_dist)
logical  :: conv(n_dist)
integer  :: rank_logl(n_dist), rank_aic(n_dist), rank_bic(n_dist)

! working scalars
real(dp) :: omega, alpha, beta, theta, h_unc, vol_ann, logl, aic, bic
real(dp) :: nu, gamma_sk, lam, mu1, mu2, sig1, sig2

character(len=*), parameter :: prices_file = "vix_spy.csv"

! ── format strings ────────────────────────────────────────────────────────────
! columns: Dist(9), omega(12), alpha(8), theta(9), beta(9),
!          d1(8), d2(8), d3(8), d4(8), mu1(8), mu2(8),
!          vol%(9), logL(12), AIC(12), BIC(12), iter(6), skew(8), ekurt(8), 3×rank(6)
! total width = 9+12+8+9+9+8+8+8+8+8+8+9+12+12+12+6+8+8+6+6+6 = 180
character(len=*), parameter :: hdr = &
    "(A9,A12,A8,A9,A9,A8,A8,A8,A8,A8,A8,A9,A12,A12,A12,A6,A8,A8,A6,A6,A6)"
character(len=*), parameter :: fmt_n  = &
    "(A9,ES12.3,F8.4,F9.4,F9.4,8X,8X,8X,8X,8X,8X,F9.2,F12.1,F12.1,F12.1,I6,F8.3,F8.3,3I6)"
character(len=*), parameter :: fmt_1p = &
    "(A9,ES12.3,F8.4,F9.4,F9.4,F8.3,8X,8X,8X,8X,8X,F9.2,F12.1,F12.1,F12.1,I6,F8.3,F8.3,3I6)"
character(len=*), parameter :: fmt_2p = &
    "(A9,ES12.3,F8.4,F9.4,F9.4,F8.3,F8.3,8X,8X,8X,8X,F9.2,F12.1,F12.1,F12.1,I6,F8.3,F8.3,3I6)"
character(len=*), parameter :: fmt_3p = &
    "(A9,ES12.3,F8.4,F9.4,F9.4,F8.3,F8.3,F8.3,8X,F8.4,F8.4,F9.2,F12.1,F12.1,F12.1,I6,F8.3,F8.3,3I6)"
character(len=*), parameter :: fmt_4p = &
    "(A9,ES12.3,F8.4,F9.4,F9.4,F8.3,F8.3,F8.3,F8.2,F8.4,F8.4,F9.2,F12.1,F12.1,F12.1,I6,F8.3,F8.3,3I6)"

call cpu_time(t_start)

call read_price_csv(prices_file, dates, col_names, prices)
nprices = size(prices, 1)
ncols   = size(col_names)
nall    = nprices - 1
nobs    = min(nret, nall)
i1      = nprices - nobs
allocate(raw_ret(nobs), ret(nobs))

write(*, '(A,I0,A,I0,A)') "Using last ", nobs, " of ", nall, " observations"
write(*, *)

do icol = 1, ncols

    raw_ret  = log(prices(i1+1:nprices, icol) / prices(i1:nprices-1, icol))
    ret_mean = mean(raw_ret)
    ret      = raw_ret - ret_mean
    ret_std  = sd(ret)

    write(*, '(A,A)') "Asset: ", trim(col_names(icol))
    write(*, '(A,F8.4,A,F8.4,A)') &
        "  mean (bps): ", ret_mean*1.0e4_dp, &
        "   emp vol_ann%: ", ret_std*sqrt(trading_days)*100.0_dp, "%"
    write(*, hdr) "Dist", "omega", "alpha", "theta", "beta", &
        "d1", "d2", "d3", "d4", "mu1", "mu2", "vol_ann%", "logL", "AIC", "BIC", &
        "iter", "skew", "ekurt", "rklogL", "rkAIC", "rkBIC"
    write(*, '(A)') repeat("-", 180)

    ! ── Normal ───────────────────────────────────────────────────────────────
    id = 1
    call flex_set_data(ret, nobs)
    call flex_set_types(proc_nagarch, dist_normal)
    np = flex_np()
    allocate(p(np), p0(np))
    call nagarch_inv_transform(1.0e-5_dp, 0.06_dp, 0.88_dp, 0.5_dp, p0)
    p = p0
    call bfgs_minimize(flex_obj, p, np, max_iter, gtol, fopt, niter, converged)
    call nagarch_transform(p(1:4), omega, alpha, beta, theta)
    h_unc   = omega / max(1.0_dp - alpha*(1.0_dp+theta**2) - beta, 1.0e-8_dp)
    vol_ann = sqrt(trading_days * h_unc) * 100.0_dp
    logl = -nobs*fopt;  aic = 2.0_dp*np - 2.0_dp*logl;  bic = np*log(real(nobs,dp)) - 2.0_dp*logl
    omegas(id) = omega;  alphas(id) = alpha;  betas(id) = beta;  thetas(id) = theta
    d1s(id) = 0.0_dp;  d2s(id) = 0.0_dp;  d3s(id) = 0.0_dp;  d4s(id) = 0.0_dp
    mu1s(id) = 0.0_dp;  mu2s(id) = 0.0_dp
    vol_anns(id) = vol_ann;  logls(id) = logl;  aics(id) = aic;  bics(id) = bic
    niters(id) = niter;  conv(id) = converged
    call flex_skew_kurt(p, np, skews(id), kurts(id))
    deallocate(p, p0)

    ! ── Student-t ─────────────────────────────────────────────────────────────
    id = 2
    call flex_set_types(proc_nagarch, dist_t)
    np = flex_np()
    allocate(p(np), p0(np))
    call nagarch_inv_transform(1.0e-5_dp, 0.06_dp, 0.88_dp, 0.5_dp, p0(1:4))
    p0(np) = p_t0
    p = p0
    call bfgs_minimize(flex_obj, p, np, max_iter, gtol, fopt, niter, converged)
    call nagarch_transform(p(1:4), omega, alpha, beta, theta)
    nu      = 2.0_dp + 98.0_dp / (1.0_dp + exp(-p(np)))
    h_unc   = omega / max(1.0_dp - alpha*(1.0_dp+theta**2) - beta, 1.0e-8_dp)
    vol_ann = sqrt(trading_days * h_unc) * 100.0_dp
    logl = -nobs*fopt;  aic = 2.0_dp*np - 2.0_dp*logl;  bic = np*log(real(nobs,dp)) - 2.0_dp*logl
    omegas(id) = omega;  alphas(id) = alpha;  betas(id) = beta;  thetas(id) = theta
    d1s(id) = nu;  d2s(id) = 0.0_dp;  d3s(id) = 0.0_dp;  d4s(id) = 0.0_dp
    mu1s(id) = 0.0_dp;  mu2s(id) = 0.0_dp
    vol_anns(id) = vol_ann;  logls(id) = logl;  aics(id) = aic;  bics(id) = bic
    niters(id) = niter;  conv(id) = converged
    call flex_skew_kurt(p, np, skews(id), kurts(id))
    deallocate(p, p0)

    ! ── Skewed-t ──────────────────────────────────────────────────────────────
    id = 3
    call flex_set_types(proc_nagarch, dist_skew_t)
    np = flex_np()
    allocate(p(np), p0(np))
    call nagarch_inv_transform(1.0e-5_dp, 0.06_dp, 0.88_dp, 0.5_dp, p0(1:4))
    p0(np-1) = p_t0;  p0(np) = 0.0_dp
    p = p0
    call bfgs_minimize(flex_obj, p, np, max_iter, gtol, fopt, niter, converged)
    call nagarch_transform(p(1:4), omega, alpha, beta, theta)
    nu       = 2.0_dp + 98.0_dp / (1.0_dp + exp(-p(np-1)))
    gamma_sk = exp(p(np))
    h_unc   = omega / max(1.0_dp - alpha*(1.0_dp+theta**2) - beta, 1.0e-8_dp)
    vol_ann = sqrt(trading_days * h_unc) * 100.0_dp
    logl = -nobs*fopt;  aic = 2.0_dp*np - 2.0_dp*logl;  bic = np*log(real(nobs,dp)) - 2.0_dp*logl
    omegas(id) = omega;  alphas(id) = alpha;  betas(id) = beta;  thetas(id) = theta
    d1s(id) = nu;  d2s(id) = gamma_sk;  d3s(id) = 0.0_dp;  d4s(id) = 0.0_dp
    mu1s(id) = 0.0_dp;  mu2s(id) = 0.0_dp
    vol_anns(id) = vol_ann;  logls(id) = logl;  aics(id) = aic;  bics(id) = bic
    niters(id) = niter;  conv(id) = converged
    call flex_skew_kurt(p, np, skews(id), kurts(id))
    deallocate(p, p0)

    ! ── NM2 ───────────────────────────────────────────────────────────────────
    id = 4
    np = nm2_nparams
    allocate(p(np), p0(np))
    call nm2_set_data(ret, nobs)
    call nagarch_inv_transform(1.0e-5_dp, 0.06_dp, 0.88_dp, 0.5_dp, p0(1:4))
    p0(5) = 0.0_dp   ! lambda = 0.5
    p0(6) = 0.0_dp   ! delta  = 0   -> mu1=mu2=0 -> N(0,1) start
    p0(7) = 0.0_dp   ! s      = 0.5 -> sig1=sig2=1
    p = p0
    call bfgs_minimize(nm2_obj, p, np, max_iter, gtol, fopt, niter, converged)
    call nagarch_transform(p(1:4), omega, alpha, beta, theta)
    call nm2_unpack_dist(p(5), p(6), p(7), lam, mu1, mu2, sig1, sig2)
    h_unc   = omega / max(1.0_dp - alpha*(1.0_dp+theta**2) - beta, 1.0e-8_dp)
    vol_ann = sqrt(trading_days * h_unc) * 100.0_dp
    logl = -nobs*fopt;  aic = 2.0_dp*np - 2.0_dp*logl;  bic = np*log(real(nobs,dp)) - 2.0_dp*logl
    omegas(id) = omega;  alphas(id) = alpha;  betas(id) = beta;  thetas(id) = theta
    d1s(id) = lam;  d2s(id) = sig1;  d3s(id) = sig2;  d4s(id) = 0.0_dp
    mu1s(id) = mu1;  mu2s(id) = mu2
    vol_anns(id) = vol_ann;  logls(id) = logl;  aics(id) = aic;  bics(id) = bic
    niters(id) = niter;  conv(id) = converged
    call nm2_skew_kurt(p, skews(id), kurts(id))
    deallocate(p, p0)

    ! ── TM2 ───────────────────────────────────────────────────────────────────
    id = 5
    np = tm2_nparams
    allocate(p(np), p0(np))
    call tm2_set_data(ret, nobs)
    call nagarch_inv_transform(1.0e-5_dp, 0.06_dp, 0.88_dp, 0.5_dp, p0(1:4))
    p0(5) = 0.0_dp   ! lambda = 0.5
    p0(6) = 0.0_dp   ! delta  = 0
    p0(7) = 0.0_dp   ! s      = 0.5
    p0(8) = p_t0     ! nu ~ 8
    p = p0
    call bfgs_minimize(tm2_obj, p, np, max_iter, gtol, fopt, niter, converged)
    call nagarch_transform(p(1:4), omega, alpha, beta, theta)
    call tm2_unpack_dist(p(5), p(6), p(7), p(8), lam, mu1, mu2, sig1, sig2, nu)
    h_unc   = omega / max(1.0_dp - alpha*(1.0_dp+theta**2) - beta, 1.0e-8_dp)
    vol_ann = sqrt(trading_days * h_unc) * 100.0_dp
    logl = -nobs*fopt;  aic = 2.0_dp*np - 2.0_dp*logl;  bic = np*log(real(nobs,dp)) - 2.0_dp*logl
    omegas(id) = omega;  alphas(id) = alpha;  betas(id) = beta;  thetas(id) = theta
    d1s(id) = lam;  d2s(id) = sig1;  d3s(id) = sig2;  d4s(id) = nu
    mu1s(id) = mu1;  mu2s(id) = mu2
    vol_anns(id) = vol_ann;  logls(id) = logl;  aics(id) = aic;  bics(id) = bic
    niters(id) = niter;  conv(id) = converged
    call tm2_skew_kurt(p, skews(id), kurts(id))
    deallocate(p, p0)

    ! ── rank and print ────────────────────────────────────────────────────────
    call rank_desc(logls, rank_logl)
    call rank_asc(aics,   rank_aic)
    call rank_asc(bics,   rank_bic)

    do id = 1, n_dist
        select case (id)
        case (1)   ! Normal: no dist params
            write(*, fmt_n) trim(dist_names(id)), &
                omegas(id), alphas(id), thetas(id), betas(id), &
                vol_anns(id), logls(id), aics(id), bics(id), niters(id), &
                skews(id), kurts(id), rank_logl(id), rank_aic(id), rank_bic(id)
        case (2)   ! t: d1=nu
            write(*, fmt_1p) trim(dist_names(id)), &
                omegas(id), alphas(id), thetas(id), betas(id), d1s(id), &
                vol_anns(id), logls(id), aics(id), bics(id), niters(id), &
                skews(id), kurts(id), rank_logl(id), rank_aic(id), rank_bic(id)
        case (3)   ! skew-t: d1=nu, d2=gamma
            write(*, fmt_2p) trim(dist_names(id)), &
                omegas(id), alphas(id), thetas(id), betas(id), d1s(id), d2s(id), &
                vol_anns(id), logls(id), aics(id), bics(id), niters(id), &
                skews(id), kurts(id), rank_logl(id), rank_aic(id), rank_bic(id)
        case (4)   ! NM2: d1=lambda, d2=sig1, d3=sig2, d4 blank, mu1, mu2
            write(*, fmt_3p) trim(dist_names(id)), &
                omegas(id), alphas(id), thetas(id), betas(id), &
                d1s(id), d2s(id), d3s(id), mu1s(id), mu2s(id), &
                vol_anns(id), logls(id), aics(id), bics(id), niters(id), &
                skews(id), kurts(id), rank_logl(id), rank_aic(id), rank_bic(id)
        case (5)   ! TM2: d1=lambda, d2=sig1, d3=sig2, d4=nu, mu1, mu2
            write(*, fmt_4p) trim(dist_names(id)), &
                omegas(id), alphas(id), thetas(id), betas(id), &
                d1s(id), d2s(id), d3s(id), d4s(id), mu1s(id), mu2s(id), &
                vol_anns(id), logls(id), aics(id), bics(id), niters(id), &
                skews(id), kurts(id), rank_logl(id), rank_aic(id), rank_bic(id)
        end select
        if (.not. conv(id)) write(*, '(4X,A)') "WARNING: did not converge"
    end do

    write(*, *)

end do

write(*, '(A)') "NAGARCH: h_t = omega + alpha*(y_{t-1} - theta*sqrt(h_{t-1}))^2 + beta*h_{t-1}"
write(*, '(A)') "vol_ann% = sqrt(252 * omega/(1 - alpha*(1+theta^2) - beta)) * 100"
write(*, '(A)') "t:      d1 = nu (degrees of freedom)"
write(*, '(A)') "Skt:    d1 = nu,    d2 = gamma (Fernandez-Steel skewness, gamma>1 => right skew)"
write(*, '(A)') "NM2:    d1 = lambda (mixing weight), d2 = sigma1, d3 = sigma2, mu1, mu2"
write(*, '(A)') "        f(z) = lambda*N(z;mu1,sig1) + (1-lambda)*N(z;mu2,sig2)"
write(*, '(A)') "        mu1 = delta*(1-lambda), mu2 = -delta*lambda  (zero-mean constraint)"
write(*, '(A)') "        V = 1 - delta^2*lambda*(1-lambda); sig1 = sqrt(s*V/lambda)"
write(*, '(A)') "        NM2 has 3 free dist params"
write(*, '(A)') "TM2:    d1 = lambda (mixing weight), d2 = sigma1, d3 = sigma2, d4 = nu, mu1, mu2"
write(*, '(A)') "        f(z) = lambda*t(z;mu1,sig1,nu) + (1-lambda)*t(z;mu2,sig2,nu)  (shared nu)"
write(*, '(A)') "        sig_k^2 = s_k*V*(nu-2)/(nu*lambda_k)  (scaled for unit t-mixture variance)"
write(*, '(A)') "        TM2 has 4 free dist params"

call cpu_time(t_end)
write(*, '(/,A,F0.3,A)') "elapsed time: ", t_end - t_start, " s"

end program xnagarch_mix_t

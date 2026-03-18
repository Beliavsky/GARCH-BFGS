! STNAGARCH × {Normal, t, skew-t, NM2, TM2}.
!
! Model: h_t = omega + [alpha_L*(1-G_t) + alpha_R*G_t]*u_t^2 + beta*h_{t-1}
!   u_t = r_{t-1} - theta*sqrt(h_{t-1}),   G_t = sigma(gamma*u_t)
!   gamma -> 0 : reduces to NAGARCH with alpha = (alpha_L+alpha_R)/2
!   alpha_L > alpha_R : stronger response to shocks left of NIC minimum
!   NIC minimum fixed at r = theta*sqrt(h_bar)  (same location as NAGARCH)
!
! Structural parameters p(1:6):
!   p1: log(omega)
!   p2,p3: softmax -> alpha_avg*(1+theta^2), beta  [stationarity]
!   p4: atanh(delta), delta = (alpha_R-alpha_L)/(alpha_L+alpha_R)
!   p5: theta  (NIC minimum location, free)
!   p6: log(gamma)  (transition speed, >0)
! p(7..np): distribution params  (nd = 0,1,2,3,4 for Normal,t,skew-t,NM2,TM2)

module stgarch_mod
    use kind_mod,       only: dp
    use math_const_mod, only: log_sqrt_2pi, pi
    implicit none
    private

    integer, parameter, public :: st_dist_normal = 1
    integer, parameter, public :: st_dist_t      = 2
    integer, parameter, public :: st_dist_skewt  = 3
    integer, parameter, public :: st_dist_nm2    = 4
    integer, parameter, public :: st_dist_tm2    = 5
    integer, parameter, public :: st_np_struct   = 6

    real(dp), allocatable, save :: st_obs(:)
    integer,               save :: st_nobs = 0
    integer,               save :: st_dist = st_dist_normal

    public :: st_set_data, st_set_dist, st_np, st_obj, &
              st_transform, st_inv_transform, st_skew_kurt, st_unpack_dist

contains

    subroutine st_set_data(y, n)
        integer,  intent(in) :: n
        real(dp), intent(in) :: y(n)
        if (allocated(st_obs)) deallocate(st_obs)
        allocate(st_obs(n))
        st_obs  = y
        st_nobs = n
    end subroutine st_set_data

    subroutine st_set_dist(dist)
        integer, intent(in) :: dist
        st_dist = dist
    end subroutine st_set_dist

    integer function st_np()
        select case (st_dist)
        case (st_dist_t);      st_np = 7
        case (st_dist_skewt);  st_np = 8
        case (st_dist_nm2);    st_np = 9
        case (st_dist_tm2);    st_np = 10
        case default;          st_np = 6   ! Normal (and fallback)
        end select
    end function st_np

    ! Decode p(1:6) -> structural params.
    ! 3-way softmax: e2=exp(p2), e3=exp(p3), S=1+e2+e3
    !   alpha_avg_raw = e2/S,  beta = e3/S
    !   alpha_avg = alpha_avg_raw / (1+theta^2)
    !   alpha_L = alpha_avg*(1-delta),  alpha_R = alpha_avg*(1+delta)
    !   delta = tanh(p4),  theta = p5,  gamma = exp(p6)
    subroutine st_transform(p, omega, alpha_l, alpha_r, beta, theta, gamma_st)
        real(dp), intent(in)  :: p(6)
        real(dp), intent(out) :: omega, alpha_l, alpha_r, beta, theta, gamma_st
        real(dp) :: e2, e3, S, alpha_avg_raw, alpha_avg, s2, delta
        omega    = exp(p(1))
        e2       = exp(p(2));  e3 = exp(p(3));  S = 1.0_dp + e2 + e3
        alpha_avg_raw = e2 / S
        beta     = e3 / S
        theta    = p(5)
        s2       = 1.0_dp + theta**2
        alpha_avg = alpha_avg_raw / s2
        delta    = tanh(p(4))
        alpha_l  = alpha_avg * (1.0_dp - delta)
        alpha_r  = alpha_avg * (1.0_dp + delta)
        gamma_st = exp(p(6))
    end subroutine st_transform

    ! Encode structural params -> unconstrained p(1:6).
    subroutine st_inv_transform(omega, alpha_l, alpha_r, beta, theta, gamma_st, p)
        real(dp), intent(in)  :: omega, alpha_l, alpha_r, beta, theta, gamma_st
        real(dp), intent(out) :: p(6)
        real(dp) :: alpha_avg, s2, A, slack, delta
        real(dp), parameter :: eps = 1.0e-10_dp
        alpha_avg = (alpha_l + alpha_r) * 0.5_dp
        s2        = 1.0_dp + theta**2
        A         = alpha_avg * s2          ! alpha_avg_raw
        slack     = max(1.0_dp - A - beta, eps)
        p(1)      = log(omega)
        p(2)      = log(A / slack)
        p(3)      = log(beta / slack)
        delta     = (alpha_r - alpha_l) / max(alpha_r + alpha_l, eps)
        delta     = max(min(delta, 1.0_dp - eps), -(1.0_dp - eps))
        p(4)      = 0.5_dp * log((1.0_dp + delta) / (1.0_dp - delta))
        p(5)      = theta
        p(6)      = log(gamma_st)
    end subroutine st_inv_transform

    ! Decode distribution display values from p_dist(1:nd).
    subroutine st_unpack_dist(p_dist, nd, d1, d2, d3, d4, mu1, mu2)
        integer,  intent(in)  :: nd
        real(dp), intent(in)  :: p_dist(nd)
        real(dp), intent(out) :: d1, d2, d3, d4, mu1, mu2
        real(dp) :: lam, nu, sig1, sig2
        d1 = 0.0_dp;  d2 = 0.0_dp;  d3 = 0.0_dp;  d4 = 0.0_dp
        mu1 = 0.0_dp;  mu2 = 0.0_dp
        select case (st_dist)
        case (st_dist_t)
            d1 = 2.0_dp + 98.0_dp / (1.0_dp + exp(-p_dist(1)))  ! nu
        case (st_dist_skewt)
            d1 = 2.0_dp + 98.0_dp / (1.0_dp + exp(-p_dist(1)))
            d2 = exp(p_dist(2))   ! gamma_sk
        case (st_dist_nm2)
            call sn_nm2_unpack(p_dist(1), p_dist(2), p_dist(3), lam, mu1, mu2, sig1, sig2)
            d1 = lam;  d2 = sig1;  d3 = sig2
        case (st_dist_tm2)
            call sn_tm2_unpack(p_dist(1), p_dist(2), p_dist(3), p_dist(4), &
                               lam, mu1, mu2, sig1, sig2, nu)
            d1 = lam;  d2 = sig1;  d3 = sig2;  d4 = nu
        end select
    end subroutine st_unpack_dist

    ! NLL and d(log f)/dz for the current distribution at z.
    ! dl_dz = d(log f)/dz  (score w.r.t. z; note: NOT d(nll)/dz)
    subroutine st_nll_score_z(z, p_dist, nll_z, dl_dz)
        real(dp), intent(in)  :: z
        real(dp), intent(in)  :: p_dist(*)
        real(dp), intent(out) :: nll_z, dl_dz
        real(dp) :: nu, gamma_sk, f_const
        select case (st_dist)
        case (st_dist_normal)
            nll_z = log_sqrt_2pi + 0.5_dp * z**2
            dl_dz = -z
        case (st_dist_t)
            nu      = 2.0_dp + 98.0_dp / (1.0_dp + exp(-p_dist(1)))
            f_const = -log_gamma(0.5_dp*(nu+1.0_dp)) + log_gamma(0.5_dp*nu) &
                      + 0.5_dp*log(pi*(nu-2.0_dp))
            nll_z   = f_const + 0.5_dp*(nu+1.0_dp)*log(1.0_dp + z**2/(nu-2.0_dp))
            dl_dz   = -(nu+1.0_dp)*z / (nu-2.0_dp + z**2)
        case (st_dist_skewt)
            nu      = 2.0_dp + 98.0_dp / (1.0_dp + exp(-p_dist(1)))
            gamma_sk = exp(p_dist(2))
            call sn_skewt_nll_dl(z, nu, gamma_sk, nll_z, dl_dz)
        case (st_dist_nm2)
            call sn_nm2_score(z, p_dist(1), p_dist(2), p_dist(3), nll_z, dl_dz)
        case (st_dist_tm2)
            call sn_tm2_score(z, p_dist(1), p_dist(2), p_dist(3), p_dist(4), nll_z, dl_dz)
        end select
    end subroutine st_nll_score_z

    ! STNAGARCH objective: NLL/n and gradient w.r.t. unconstrained p(1:np).
    ! Structural gradient (p1..p6): analytic via coupled dh recursion.
    ! Distribution gradient (p7..np): central differences.
    subroutine st_obj(p, np, f, g)
        integer,  intent(in)  :: np
        real(dp), intent(in)  :: p(np)
        real(dp), intent(out) :: f, g(np)

        real(dp) :: omega, alpha_l, alpha_r, beta, theta, gamma_st
        real(dp) :: s2, A, B, D, alpha_avg
        real(dp) :: h, sqrth, u, Gval, gdot, alpha_t
        real(dp) :: z, nll_z, dl_dz, factor, kappa
        real(dp) :: dh(6)          ! d(h_t)/d(omega,alpha_L,alpha_R,beta,theta,gamma)
        real(dp) :: dir(6)         ! direct contribution at current t
        real(dp) :: gs(6)          ! structural gradient accumulators
        real(dp) :: p_dist_tmp(10) ! workspace for FD perturbation (max 4 dist params)
        real(dp) :: dummy_p(1)     ! placeholder for Normal (no dist params)
        real(dp) :: nll_p, nll_m, dl_dummy
        real(dp) :: fac_alpha, delta
        real(dp), parameter :: h_fd  = 1.0e-5_dp
        real(dp), parameter :: eps_h = 1.0e-12_dp
        integer :: t, nd, k

        nd = np - st_np_struct

        call st_transform(p(1:6), omega, alpha_l, alpha_r, beta, theta, gamma_st)
        s2       = 1.0_dp + theta**2
        alpha_avg = (alpha_l + alpha_r) * 0.5_dp
        A        = alpha_avg * s2            ! alpha_avg_raw = alpha_avg*(1+theta^2)
        B        = beta
        D        = max(1.0_dp - A - B, 1.0e-8_dp)

        ! Unconditional variance as h_1
        h        = omega / D

        ! Initial dh/d(structural) from stationarity approximation
        dh(1) = 1.0_dp / D
        dh(2) = h * s2 / (2.0_dp * D)         ! d h_bar / d alpha_L
        dh(3) = h * s2 / (2.0_dp * D)         ! d h_bar / d alpha_R
        dh(4) = h / D                          ! d h_bar / d beta
        dh(5) = h * alpha_avg * 2.0_dp*theta / D  ! d h_bar / d theta
        dh(6) = 0.0_dp                         ! d h_bar / d gamma ~ 0

        f  = 0.0_dp
        gs = 0.0_dp
        g(st_np_struct+1:np) = 0.0_dp
        dummy_p(1) = 0.0_dp

        do t = 1, st_nobs
            sqrth = sqrt(max(h, eps_h))
            z     = st_obs(t) / sqrth

            ! NLL and score; for Normal nd=0 so p_dist is never read
            if (nd > 0) then
                call st_nll_score_z(z, p(st_np_struct+1), nll_z, dl_dz)
            else
                call st_nll_score_z(z, dummy_p, nll_z, dl_dz)
            end if
            factor = (1.0_dp + dl_dz * z) / (2.0_dp * max(h, eps_h))
            f      = f + nll_z + 0.5_dp * log(max(h, eps_h))
            gs(:)  = gs(:) + factor * dh(:)

            ! Central-difference gradient for distribution params
            if (nd > 0) then
                p_dist_tmp(1:nd) = p(st_np_struct+1:np)
                do k = 1, nd
                    p_dist_tmp(k) = p_dist_tmp(k) + h_fd
                    call st_nll_score_z(z, p_dist_tmp(1), nll_p, dl_dummy)
                    p_dist_tmp(k) = p_dist_tmp(k) - 2.0_dp*h_fd
                    call st_nll_score_z(z, p_dist_tmp(1), nll_m, dl_dummy)
                    g(st_np_struct+k) = g(st_np_struct+k) + (nll_p - nll_m) / (2.0_dp * h_fd)
                    p_dist_tmp(k) = p_dist_tmp(k) + h_fd   ! restore
                end do
            end if

            ! Advance STNAGARCH recurrence
            u       = st_obs(t) - theta * sqrth
            Gval    = 1.0_dp / (1.0_dp + exp(-gamma_st * u))
            gdot    = Gval * (1.0_dp - Gval)
            alpha_t = alpha_l*(1.0_dp - Gval) + alpha_r*Gval

            ! Direct contributions dir(k) = d(h_{t+1})/d(param_k) holding h_t fixed
            dir(1) = 1.0_dp
            dir(2) = (1.0_dp - Gval) * u**2
            dir(3) = Gval * u**2
            dir(4) = h
            dir(5) = -sqrth * (2.0_dp*alpha_t*u + (alpha_r-alpha_l)*gdot*gamma_st*u**2)
            dir(6) = (alpha_r - alpha_l) * gdot * u**3

            ! kappa = d(h_{t+1})/d(h_t)
            kappa  = beta - (theta/sqrth) * &
                     (alpha_t*u + (alpha_r-alpha_l)*gdot*gamma_st*0.5_dp*u**2)

            dh(:) = dir(:) + kappa * dh(:)
            h     = omega + alpha_t * u**2 + beta * h
        end do

        f  = f  / st_nobs
        gs = gs / st_nobs
        g(st_np_struct+1:np) = g(st_np_struct+1:np) / st_nobs

        ! Chain rule: structural (omega,alpha_L,alpha_R,beta,theta,gamma) -> p(1:6)
        ! Decode softmax state
        delta     = tanh(p(4))
        A         = alpha_avg * s2    ! e2/S
        B         = beta              ! e3/S (= beta)
        fac_alpha = (gs(2)*(1.0_dp-delta) + gs(3)*(1.0_dp+delta)) / s2

        g(1) = gs(1) * omega
        g(2) = fac_alpha * A*(1.0_dp-A) - gs(4)*A*B               ! softmax A: dA/dp2 = A*(1-A)
        g(3) = fac_alpha * (-A*B) + gs(4)*B*(1.0_dp-B)           ! softmax B: dB/dp3 = B*(1-B)
        g(4) = (gs(3) - gs(2)) * alpha_avg * (1.0_dp - delta**2)  ! tanh
        g(5) = -2.0_dp*theta/s2 * (gs(2)*alpha_l + gs(3)*alpha_r) + gs(5)  ! theta
        g(6) = gs(6) * gamma_st
    end subroutine st_obj

    ! Empirical skewness and excess kurtosis of standardised residuals z_t = y_t/sqrt(h_t).
    subroutine st_skew_kurt(p, skew, kurt)
        real(dp), intent(in)  :: p(*)
        real(dp), intent(out) :: skew, kurt
        real(dp) :: omega, alpha_l, alpha_r, beta, theta, gamma_st
        real(dp) :: h, sqrth, u, Gval, alpha_t
        real(dp), allocatable :: zz(:)
        real(dp) :: zm, zv, zs, zk, dz, dz2, rn
        integer  :: t
        call st_transform(p(1:6), omega, alpha_l, alpha_r, beta, theta, gamma_st)
        h = omega / max(1.0_dp - (alpha_l+alpha_r)*0.5_dp*(1.0_dp+theta**2) - beta, 1.0e-8_dp)
        allocate(zz(st_nobs))
        do t = 1, st_nobs
            sqrth  = sqrt(max(h, 1.0e-12_dp))
            zz(t)  = st_obs(t) / sqrth
            u       = st_obs(t) - theta * sqrth
            Gval    = 1.0_dp / (1.0_dp + exp(-gamma_st * u))
            alpha_t = alpha_l*(1.0_dp-Gval) + alpha_r*Gval
            h      = omega + alpha_t * u**2 + beta * h
        end do
        rn = real(st_nobs, dp)
        zm = sum(zz) / rn
        zv = 0.0_dp;  zs = 0.0_dp;  zk = 0.0_dp
        do t = 1, st_nobs
            dz = zz(t) - zm;  dz2 = dz**2
            zv = zv + dz2;  zs = zs + dz*dz2;  zk = zk + dz2**2
        end do
        zv   = zv / rn
        skew = (zs/rn) / zv**1.5_dp
        kurt = (zk/rn) / zv**2 - 3.0_dp
        deallocate(zz)
    end subroutine st_skew_kurt

    ! ── Fernandez-Steel skewed-t (same as in garch_flex_mod) ─────────────────
    subroutine sn_skewt_nll_dl(z, nu, gamma, nll, dl_dz)
        real(dp), intent(in)  :: z, nu, gamma
        real(dp), intent(out) :: nll, dl_dz
        real(dp) :: c, mu1, sig2, sig, x, u, du_dz, log1pu2
        c    = exp(log_gamma(0.5_dp*(nu+1.0_dp)) - log_gamma(0.5_dp*nu)) &
               / sqrt(pi*(nu-2.0_dp))
        mu1  = 2.0_dp*(gamma - 1.0_dp/gamma)*c*(nu-2.0_dp)/(nu-1.0_dp)
        sig2 = gamma**2 + 1.0_dp/gamma**2 - 1.0_dp - mu1**2
        sig  = sqrt(sig2)
        x    = sig*z + mu1
        if (x >= 0.0_dp) then
            u     = x / gamma
            du_dz = sig / gamma
        else
            u     = x * gamma
            du_dz = sig * gamma
        end if
        log1pu2 = log(1.0_dp + u**2/(nu-2.0_dp))
        nll = log(0.5_dp*(gamma + 1.0_dp/gamma)) &
              + log_gamma(0.5_dp*nu) - log_gamma(0.5_dp*(nu+1.0_dp)) &
              + 0.5_dp*log(pi*(nu-2.0_dp)) &
              + 0.5_dp*(nu+1.0_dp)*log1pu2 &
              + log(sig)
        dl_dz = -(nu+1.0_dp)*u/(nu-2.0_dp + u**2) * du_dz
    end subroutine sn_skewt_nll_dl

    ! ── NM2 helpers ──────────────────────────────────────────────────────────
    subroutine sn_nm2_unpack(p5, p6, p7, lam, mu1, mu2, sig1, sig2)
        real(dp), intent(in)  :: p5, p6, p7
        real(dp), intent(out) :: lam, mu1, mu2, sig1, sig2
        real(dp) :: delta, s, V
        real(dp), parameter :: eps = 1.0e-8_dp
        lam   = 1.0_dp / (1.0_dp + exp(-p5))
        delta = p6;  s = 1.0_dp / (1.0_dp + exp(-p7))
        mu1   = delta*(1.0_dp - lam);  mu2 = -delta*lam
        V     = max(1.0_dp - delta**2*lam*(1.0_dp-lam), eps)
        sig1  = sqrt(s*V / max(lam, eps))
        sig2  = sqrt((1.0_dp-s)*V / max(1.0_dp-lam, eps))
    end subroutine sn_nm2_unpack

    subroutine sn_nm2_score(z, p5, p6, p7, nll_z, dl_dz)
        real(dp), intent(in)  :: z, p5, p6, p7
        real(dp), intent(out) :: nll_z, dl_dz
        real(dp) :: lam, mu1, mu2, sig1, sig2, u1, u2, phi1, phi2, w1, w2, f_mix
        real(dp), parameter :: eps = 1.0e-300_dp
        call sn_nm2_unpack(p5, p6, p7, lam, mu1, mu2, sig1, sig2)
        u1 = (z-mu1)/sig1;  u2 = (z-mu2)/sig2
        phi1 = exp(-0.5_dp*u1**2)/sig1;  phi2 = exp(-0.5_dp*u2**2)/sig2
        w1 = lam*phi1;  w2 = (1.0_dp-lam)*phi2
        f_mix = max(w1+w2, eps)
        nll_z = log_sqrt_2pi - log(f_mix)
        dl_dz = -(w1*u1/sig1 + w2*u2/sig2) / f_mix
    end subroutine sn_nm2_score

    ! ── TM2 helpers ──────────────────────────────────────────────────────────
    subroutine sn_tm2_unpack(p5, p6, p7, p8, lam, mu1, mu2, sig1, sig2, nu)
        real(dp), intent(in)  :: p5, p6, p7, p8
        real(dp), intent(out) :: lam, mu1, mu2, sig1, sig2, nu
        real(dp) :: delta, s, V, scale
        real(dp), parameter :: eps = 1.0e-8_dp
        lam   = 1.0_dp / (1.0_dp + exp(-p5))
        delta = p6;  s = 1.0_dp / (1.0_dp + exp(-p7))
        nu    = 2.0_dp + 98.0_dp / (1.0_dp + exp(-p8))
        mu1   = delta*(1.0_dp-lam);  mu2 = -delta*lam
        V     = max(1.0_dp - delta**2*lam*(1.0_dp-lam), eps)
        scale = (nu-2.0_dp)/nu
        sig1  = sqrt(s*V*scale / max(lam, eps))
        sig2  = sqrt((1.0_dp-s)*V*scale / max(1.0_dp-lam, eps))
    end subroutine sn_tm2_unpack

    subroutine sn_tm2_score(z, p5, p6, p7, p8, nll_z, dl_dz)
        real(dp), intent(in)  :: z, p5, p6, p7, p8
        real(dp), intent(out) :: nll_z, dl_dz
        real(dp) :: lam, mu1, mu2, sig1, sig2, nu, u1, u2, g1, g2, w1, w2, f_mix, log_ct
        real(dp), parameter :: eps = 1.0e-300_dp
        call sn_tm2_unpack(p5, p6, p7, p8, lam, mu1, mu2, sig1, sig2, nu)
        log_ct = log_gamma(0.5_dp*(nu+1.0_dp)) - log_gamma(0.5_dp*nu) &
                 - 0.5_dp*log(acos(-1.0_dp)*nu)
        u1 = (z-mu1)/sig1;  u2 = (z-mu2)/sig2
        g1 = (1.0_dp + u1**2/nu)**(-(nu+1.0_dp)*0.5_dp) / sig1
        g2 = (1.0_dp + u2**2/nu)**(-(nu+1.0_dp)*0.5_dp) / sig2
        w1 = lam*g1;  w2 = (1.0_dp-lam)*g2
        f_mix = max(w1+w2, eps)
        nll_z = -(log_ct + log(f_mix))
        dl_dz = -(w1*(nu+1.0_dp)*u1/sig1/(nu+u1**2) + w2*(nu+1.0_dp)*u2/sig2/(nu+u2**2)) / f_mix
    end subroutine sn_tm2_score

end module stgarch_mod


! ── Main program ─────────────────────────────────────────────────────────────
program xstgarch
! STNAGARCH with Normal, Student-t, skewed-t, NM2, TM2 noise distributions.

use kind_mod,      only: dp
use csv_mod,       only: read_price_csv
use stgarch_mod,   only: st_set_data, st_set_dist, st_np, st_obj, &
                          st_transform, st_inv_transform, st_skew_kurt, st_unpack_dist, &
                          st_dist_normal, st_dist_t, st_dist_skewt, st_dist_nm2, st_dist_tm2
use garch_flex_mod, only: flex_set_data, flex_set_types, flex_np, flex_obj, &
                          proc_nagarch, &
                          flex_dist_normal => dist_normal, &
                          flex_dist_t      => dist_t
use nagarch_module, only: nagarch_transform, nagarch_inv_transform
use bfgs_module,   only: bfgs_minimize
use stats_mod,     only: mean, sd
use rank_mod,      only: rank_desc, rank_asc
implicit none

integer,  parameter :: n_dist       = 5
integer,  parameter :: n_nagarch    = 2               ! NAGARCH comparison: Normal, t
integer,  parameter :: n_comp       = n_nagarch + n_dist
real(dp), parameter :: trading_days = 252.0_dp
integer,  parameter :: max_iter     = 500
real(dp), parameter :: gtol         = 1.0e-6_dp
integer,  parameter :: nret         = 10**6
real(dp), parameter :: p_t0         = -2.729_dp   ! starting logit for nu ~ 8

! Starting structural values (symmetric: alpha_L = alpha_R, delta = 0, gamma = 1)
real(dp), parameter :: st_omega0  = 1.0e-5_dp
real(dp), parameter :: st_alpha0  = 0.06_dp       ! alpha_L = alpha_R
real(dp), parameter :: st_beta0   = 0.88_dp
real(dp), parameter :: st_theta0  = 0.5_dp
real(dp), parameter :: st_gamma0  = 1.0_dp

character(len=9), parameter :: dist_names(n_dist) = &
    ["Normal   ", "Student-t", "Skewed-t ", "NM2      ", "TM2      "]

integer, parameter :: dist_ids(n_dist) = &
    [st_dist_normal, st_dist_t, st_dist_skewt, st_dist_nm2, st_dist_tm2]

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
integer  :: np, niter, id, ic
logical  :: converged
real(dp) :: t_start, t_end

! per-model results
real(dp) :: omegas(n_dist), alpha_ls(n_dist), alpha_rs(n_dist)
real(dp) :: betas(n_dist), thetas(n_dist), gammas(n_dist)
real(dp) :: d1s(n_dist), d2s(n_dist), d3s(n_dist), d4s(n_dist)
real(dp) :: mu1s(n_dist), mu2s(n_dist)
real(dp) :: vol_anns(n_dist), logls(n_dist), aics(n_dist), bics(n_dist)
real(dp) :: skews(n_dist), kurts(n_dist)
integer  :: niters(n_dist)
logical  :: conv(n_dist)
integer  :: rank_logl(n_dist), rank_aic(n_dist), rank_bic(n_dist)

! working scalars
real(dp) :: omega, alpha_l, alpha_r, beta, theta, gamma_st, h_unc, vol_ann, logl, aic, bic
real(dp) :: alpha_avg, nu_val
real(dp) :: p0_struct(6)

! compact comparison table (NAGARCH × 2 + STNAGARCH × 5)
real(dp) :: cmp_omega(n_comp), cmp_alpha(n_comp), cmp_theta(n_comp)
real(dp) :: cmp_beta(n_comp),  cmp_gamma(n_comp), cmp_shape(n_comp)
real(dp) :: cmp_vol(n_comp),   cmp_logl(n_comp),  cmp_aic(n_comp), cmp_bic(n_comp)
integer  :: cmp_niter(n_comp), cmp_rl(n_comp), cmp_ra(n_comp), cmp_rb(n_comp)
logical  :: cmp_conv(n_comp)
character(len=9) :: cmp_proc(n_comp), cmp_dist_nm(n_comp)

character(len=*), parameter :: prices_file = "vix_spy.csv"

! ── format strings ──────────────────────────────────────────────────────────
! columns: Dist(9), omega(12), alpha_L(8), alpha_R(8), theta(9), beta(9), gamma(8),
!          d1(8), d2(8), d3(8), d4(8), mu1(8), mu2(8),
!          vol%(9), logL(12), AIC(12), BIC(12), iter(6), skew(8), ekurt(8), 3xrank(6)
! total: 9+12+8+8+9+9+8+8+8+8+8+8+8+9+12+12+12+6+8+8+6+6+6 = 196
character(len=*), parameter :: hdr = &
    "(A9,A12,A8,A8,A9,A9,A8,A8,A8,A8,A8,A8,A8,A9,A12,A12,A12,A6,A8,A8,A6,A6,A6)"
character(len=*), parameter :: fmt_n  = &
    "(A9,ES12.3,F8.4,F8.4,F9.4,F9.4,F8.3,8X,8X,8X,8X,8X,8X,F9.2,F12.1,F12.1,F12.1,I6,F8.3,F8.3,3I6)"
character(len=*), parameter :: fmt_1p = &
    "(A9,ES12.3,F8.4,F8.4,F9.4,F9.4,F8.3,F8.3,8X,8X,8X,8X,8X,F9.2,F12.1,F12.1,F12.1,I6,F8.3,F8.3,3I6)"
character(len=*), parameter :: fmt_2p = &
    "(A9,ES12.3,F8.4,F8.4,F9.4,F9.4,F8.3,F8.3,F8.3,8X,8X,8X,8X,F9.2,F12.1,F12.1,F12.1,I6,F8.3,F8.3,3I6)"
character(len=*), parameter :: fmt_3p = &
    "(A9,ES12.3,F8.4,F8.4,F9.4,F9.4,F8.3,F8.3,F8.3,F8.3,8X,F8.4,F8.4,F9.2,F12.1,F12.1,F12.1,I6,F8.3,F8.3,3I6)"
character(len=*), parameter :: fmt_4p = &
    "(A9,ES12.3,F8.4,F8.4,F9.4,F9.4,F8.3,F8.3,F8.3,F8.3,F8.2,F8.4,F8.4,F9.2,F12.1,F12.1,F12.1,I6,F8.3,F8.3,3I6)"

! compact comparison table format: Process(9),1X,Dist(9),omega(12),alpha(10),theta(9),beta(9),gamma(8),shape(7),vol%(9),logL(12),AIC(12),BIC(12),iter(6),rklogL(7),rkAIC(6),rkBIC(6) = 144
character(len=*), parameter :: hdr_cmp = &
    "(A9,1X,A9,A12,A10,A9,A9,A8,A7,A9,A12,A12,A12,A6,A7,A6,A6)"
character(len=*), parameter :: fmt_cmp = &
    "(A9,1X,A9,ES12.3,F10.4,F9.4,F9.4,F8.3,F7.2,F9.2,F12.1,F12.1,F12.1,I6,I7,I6,I6)"

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

! Compute starting structural p0_struct (same for all distributions)
call st_inv_transform(st_omega0, st_alpha0, st_alpha0, st_beta0, st_theta0, st_gamma0, p0_struct)

do icol = 1, ncols

    raw_ret  = log(prices(i1+1:nprices, icol) / prices(i1:nprices-1, icol))
    ret_mean = mean(raw_ret)
    ret      = raw_ret - ret_mean
    ret_std  = sd(ret)

    write(*, '(A,A)') "Asset: ", trim(col_names(icol))
    write(*, '(A,F8.4,A,F8.4,A)') &
        "  mean (bps): ", ret_mean*1.0e4_dp, &
        "   emp vol_ann%: ", ret_std*sqrt(trading_days)*100.0_dp, "%"
    write(*, hdr) "Dist", "omega", "alpha_L", "alpha_R", "theta", "beta", "gamma", &
        "d1", "d2", "d3", "d4", "mu1", "mu2", "vol_ann%", "logL", "AIC", "BIC", &
        "iter", "skew", "ekurt", "rklogL", "rkAIC", "rkBIC"
    write(*, '(A)') repeat("-", 196)

    call st_set_data(ret, nobs)

    do id = 1, n_dist

        call st_set_dist(dist_ids(id))
        np = st_np()
        allocate(p(np), p0(np))

        p0(1:6) = p0_struct
        select case (dist_ids(id))
        case (st_dist_normal)
            ! no extra params
        case (st_dist_t)
            p0(7) = p_t0
        case (st_dist_skewt)
            p0(7) = p_t0;  p0(8) = 0.0_dp
        case (st_dist_nm2)
            p0(7) = 0.0_dp;  p0(8) = 0.0_dp;  p0(9) = 0.0_dp
        case (st_dist_tm2)
            p0(7) = 0.0_dp;  p0(8) = 0.0_dp;  p0(9) = 0.0_dp;  p0(10) = p_t0
        end select
        p = p0

        call bfgs_minimize(st_obj, p, np, max_iter, gtol, fopt, niter, converged)
        call st_transform(p(1:6), omega, alpha_l, alpha_r, beta, theta, gamma_st)
        call st_unpack_dist(p(7:np), np-6, d1s(id), d2s(id), d3s(id), d4s(id), mu1s(id), mu2s(id))

        h_unc   = omega / max(1.0_dp - (alpha_l+alpha_r)*0.5_dp*(1.0_dp+theta**2) - beta, 1.0e-8_dp)
        vol_ann = sqrt(trading_days * h_unc) * 100.0_dp
        logl    = -nobs * fopt
        aic     = 2.0_dp*np             - 2.0_dp*logl
        bic     = np*log(real(nobs,dp)) - 2.0_dp*logl

        omegas(id)   = omega;   alpha_ls(id) = alpha_l;  alpha_rs(id) = alpha_r
        betas(id)    = beta;    thetas(id)   = theta;    gammas(id)   = gamma_st
        vol_anns(id) = vol_ann; logls(id)    = logl
        aics(id)     = aic;     bics(id)     = bic
        niters(id)   = niter;   conv(id)     = converged
        call st_skew_kurt(p, skews(id), kurts(id))
        deallocate(p, p0)

    end do

    call rank_desc(logls, rank_logl)
    call rank_asc(aics,   rank_aic)
    call rank_asc(bics,   rank_bic)

    do id = 1, n_dist
        select case (dist_ids(id))
        case (st_dist_normal)
            write(*, fmt_n) trim(dist_names(id)), &
                omegas(id), alpha_ls(id), alpha_rs(id), thetas(id), betas(id), gammas(id), &
                vol_anns(id), logls(id), aics(id), bics(id), niters(id), &
                skews(id), kurts(id), rank_logl(id), rank_aic(id), rank_bic(id)
        case (st_dist_t)
            write(*, fmt_1p) trim(dist_names(id)), &
                omegas(id), alpha_ls(id), alpha_rs(id), thetas(id), betas(id), gammas(id), &
                d1s(id), &
                vol_anns(id), logls(id), aics(id), bics(id), niters(id), &
                skews(id), kurts(id), rank_logl(id), rank_aic(id), rank_bic(id)
        case (st_dist_skewt)
            write(*, fmt_2p) trim(dist_names(id)), &
                omegas(id), alpha_ls(id), alpha_rs(id), thetas(id), betas(id), gammas(id), &
                d1s(id), d2s(id), &
                vol_anns(id), logls(id), aics(id), bics(id), niters(id), &
                skews(id), kurts(id), rank_logl(id), rank_aic(id), rank_bic(id)
        case (st_dist_nm2)
            write(*, fmt_3p) trim(dist_names(id)), &
                omegas(id), alpha_ls(id), alpha_rs(id), thetas(id), betas(id), gammas(id), &
                d1s(id), d2s(id), d3s(id), mu1s(id), mu2s(id), &
                vol_anns(id), logls(id), aics(id), bics(id), niters(id), &
                skews(id), kurts(id), rank_logl(id), rank_aic(id), rank_bic(id)
        case (st_dist_tm2)
            write(*, fmt_4p) trim(dist_names(id)), &
                omegas(id), alpha_ls(id), alpha_rs(id), thetas(id), betas(id), gammas(id), &
                d1s(id), d2s(id), d3s(id), d4s(id), mu1s(id), mu2s(id), &
                vol_anns(id), logls(id), aics(id), bics(id), niters(id), &
                skews(id), kurts(id), rank_logl(id), rank_aic(id), rank_bic(id)
        end select
        if (.not. conv(id)) write(*, '(4X,A)') "WARNING: did not converge"
    end do

    ! ── NAGARCH comparison fits ───────────────────────────────────────────────
    call flex_set_data(ret, nobs)

    ! NAGARCH-Normal (np=4)
    call flex_set_types(proc_nagarch, flex_dist_normal)
    np = flex_np()
    allocate(p(np), p0(np))
    call nagarch_inv_transform(st_omega0, st_alpha0, st_beta0, st_theta0, p0(1:4))
    p = p0
    call bfgs_minimize(flex_obj, p, np, max_iter, gtol, fopt, niter, converged)
    call nagarch_transform(p(1:4), omega, alpha_l, beta, theta)
    h_unc = omega / max(1.0_dp - alpha_l*(1.0_dp+theta**2) - beta, 1.0e-8_dp)
    cmp_proc(1) = "NAGARCH  ";  cmp_dist_nm(1) = "Normal   "
    cmp_omega(1) = omega;  cmp_alpha(1) = alpha_l;  cmp_theta(1) = theta
    cmp_beta(1) = beta;    cmp_gamma(1) = 0.0_dp;   cmp_shape(1) = 0.0_dp
    cmp_vol(1)  = sqrt(trading_days * h_unc) * 100.0_dp
    cmp_logl(1) = -nobs * fopt
    cmp_aic(1)  = 2.0_dp*np - 2.0_dp*cmp_logl(1)
    cmp_bic(1)  = np*log(real(nobs,dp)) - 2.0_dp*cmp_logl(1)
    cmp_niter(1) = niter;  cmp_conv(1) = converged
    deallocate(p, p0)

    ! NAGARCH-t (np=5)
    call flex_set_types(proc_nagarch, flex_dist_t)
    np = flex_np()
    allocate(p(np), p0(np))
    call nagarch_inv_transform(st_omega0, st_alpha0, st_beta0, st_theta0, p0(1:4))
    p0(np) = p_t0
    p = p0
    call bfgs_minimize(flex_obj, p, np, max_iter, gtol, fopt, niter, converged)
    call nagarch_transform(p(1:4), omega, alpha_l, beta, theta)
    nu_val = 2.0_dp + 98.0_dp / (1.0_dp + exp(-p(np)))
    h_unc = omega / max(1.0_dp - alpha_l*(1.0_dp+theta**2) - beta, 1.0e-8_dp)
    cmp_proc(2) = "NAGARCH  ";  cmp_dist_nm(2) = "t        "
    cmp_omega(2) = omega;  cmp_alpha(2) = alpha_l;  cmp_theta(2) = theta
    cmp_beta(2) = beta;    cmp_gamma(2) = 0.0_dp;   cmp_shape(2) = nu_val
    cmp_vol(2)  = sqrt(trading_days * h_unc) * 100.0_dp
    cmp_logl(2) = -nobs * fopt
    cmp_aic(2)  = 2.0_dp*np - 2.0_dp*cmp_logl(2)
    cmp_bic(2)  = np*log(real(nobs,dp)) - 2.0_dp*cmp_logl(2)
    cmp_niter(2) = niter;  cmp_conv(2) = converged
    deallocate(p, p0)

    ! Fill STNAGARCH entries from stored section-1 results
    do id = 1, n_dist
        ic = n_nagarch + id
        cmp_proc(ic)  = "STNAGARCH"
        cmp_dist_nm(ic) = dist_names(id)
        alpha_avg     = (alpha_ls(id) + alpha_rs(id)) * 0.5_dp
        cmp_omega(ic) = omegas(id);   cmp_alpha(ic) = alpha_avg
        cmp_theta(ic) = thetas(id);   cmp_beta(ic)  = betas(id)
        cmp_gamma(ic) = gammas(id)
        select case (dist_ids(id))
        case (st_dist_t, st_dist_skewt);  cmp_shape(ic) = d1s(id)   ! nu
        case (st_dist_nm2, st_dist_tm2);  cmp_shape(ic) = d1s(id)   ! lambda
        case default;                     cmp_shape(ic) = 0.0_dp
        end select
        cmp_vol(ic)   = vol_anns(id);  cmp_logl(ic) = logls(id)
        cmp_aic(ic)   = aics(id);      cmp_bic(ic)  = bics(id)
        cmp_niter(ic) = niters(id);    cmp_conv(ic) = conv(id)
    end do

    call rank_desc(cmp_logl, cmp_rl)
    call rank_asc(cmp_aic,   cmp_ra)
    call rank_asc(cmp_bic,   cmp_rb)

    write(*, hdr_cmp) "Process", "Dist", "omega", "alpha_avg", "theta", "beta", &
        "gamma", "shape", "vol_ann%", "logL", "AIC", "BIC", "iter", "rklogL", "rkAIC", "rkBIC"
    write(*, '(A)') repeat("-", 144)
    do ic = 1, n_comp
        write(*, fmt_cmp) trim(cmp_proc(ic)), trim(cmp_dist_nm(ic)), &
            cmp_omega(ic), cmp_alpha(ic), cmp_theta(ic), cmp_beta(ic), &
            cmp_gamma(ic), cmp_shape(ic), cmp_vol(ic), &
            cmp_logl(ic), cmp_aic(ic), cmp_bic(ic), cmp_niter(ic), &
            cmp_rl(ic), cmp_ra(ic), cmp_rb(ic)
        if (.not. cmp_conv(ic)) write(*, '(4X,A)') "WARNING: did not converge"
    end do

    write(*, *)

end do

write(*, '(A)') "STNAGARCH: h_t = omega + [alpha_L*(1-G_t)+alpha_R*G_t]*u_t^2 + beta*h_{t-1}"
write(*, '(A)') "           u_t = r_{t-1} - theta*sqrt(h_{t-1}),   G_t = sigma(gamma*u_t)"
write(*, '(A)') "           NIC minimum at r = theta*sqrt(h_bar);  gamma->0: NAGARCH"
write(*, '(A)') "vol_ann% = sqrt(252 * omega/(1 - alpha_avg*(1+theta^2) - beta)) * 100"
write(*, '(A)') "t:      d1 = nu"
write(*, '(A)') "Skt:    d1 = nu,    d2 = gamma_sk (Fernandez-Steel)"
write(*, '(A)') "NM2:    d1 = lambda, d2 = sigma1, d3 = sigma2, mu1, mu2"
write(*, '(A)') "TM2:    d1 = lambda, d2 = sigma1, d3 = sigma2, d4 = nu, mu1, mu2"

call cpu_time(t_end)
write(*, '(/,A,F0.3,A)') "(2) finished xstgarch.f90, elapsed time: ", t_end - t_start, " s"

end program xstgarch

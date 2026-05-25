! GARCH(1,1) with symmetric Normal Inverse Gaussian (NIG) innovations.
!
! Innovation model:  y_t = sqrt(h_t) * epsilon_t
!   epsilon_t ~ symmetric NIG(alp): NIG(alpha=alp, beta=0, mu=0, delta=alp)
!   PDF: f(x) = (alp^2/pi) * exp(alp^2) * K1(alp*sqrt(alp^2+x^2)) / sqrt(alp^2+x^2)
!   Mean=0, Variance = delta/alpha = alp/alp = 1 for all alp > 0.
!   alp -> inf: NIG -> Normal;  small alp: heavy tails.
!
! Negative log-likelihood per observation:
!   nll_t = 0.5*log(h_t) - 2*log(alp) + log(pi) - alp^2
!           + 0.5*log(alp^2 + u_t^2) - log_k1(r_t)
!   u_t = y_t / sqrt(h_t),   r_t = alp * sqrt(alp^2 + u_t^2)
!
! Gradient w.r.t. h_t:
!   d(nll_t)/dh_t = [1 - u_t^2*(2 + r_t*K0(r_t)/K1(r_t))/(alp^2+u_t^2)] / (2*h_t)
!
! Gradient w.r.t. alp (summed over n obs, from direct differentiation):
!   d(NLL)/d(alp) = n*(-1/alp - 2*alp)
!                 + sum_t [ K0(r_t)/K1(r_t) * (2*alp^2+u_t^2)/sqrt(alp^2+u_t^2)
!                           + 3*alp/(alp^2+u_t^2) ]
!
! Unconstrained parameterisation:
!   p(1:3) -> (omega, alpha_g, beta_g) via softmax (same as Normal GARCH)
!   p(4)   -> alp = 0.1 + 19.9 * sigmoid(p4)   (alp in (0.1, 20))
!   Inverse: p4 = log((alp - 0.1) / (20 - alp))
!   d(alp)/d(p4) = (alp - 0.1) * (20 - alp) / 19.9
!
! Simulation via normal-variance mixture:
!   V ~ InvGauss(mu=1, lambda=alp^2),  Z ~ N(0,1) independent
!   epsilon = sqrt(V) * Z  (has the symmetric NIG distribution)

module garch_nig_mod
    use kind_mod,       only: dp
    use math_const_mod, only: pi
    use garch_mod,   only: garch_transform, garch_inv_transform
    use random_mod,     only: random_nig_sym
    use special_mod,    only: bessel_k01
    implicit none
    private
    real(dp), parameter :: alp_min = 0.1_dp
    real(dp), parameter :: alp_max = 20.0_dp
    real(dp), parameter :: alp_rng = alp_max - alp_min   ! 19.9

    real(dp), allocatable, save :: gni_obs(:)
    integer,               save :: gni_nobs = 0

    public :: garch_nig_set_data, garch_nig_simulate, garch_nig_obj, &
              garch_nig_transform, garch_nig_inv_transform

contains

    subroutine garch_nig_set_data(y, n)
        integer,  intent(in) :: n
        real(dp), intent(in) :: y(n)
        if (allocated(gni_obs)) deallocate(gni_obs)
        allocate(gni_obs(n))
        gni_obs  = y
        gni_nobs = n
    end subroutine garch_nig_set_data

    ! Unconstrained p(4) -> constrained (omega, alpha_g, beta_g, alp_nig).
    subroutine garch_nig_transform(p, omega, alpha_g, beta_g, alp)
        real(dp), intent(in)  :: p(4)
        real(dp), intent(out) :: omega, alpha_g, beta_g, alp
        real(dp) :: sig
        call garch_transform(p(1:3), omega, alpha_g, beta_g)
        sig = 1.0_dp / (1.0_dp + exp(-p(4)))
        alp = alp_min + alp_rng * sig
    end subroutine garch_nig_transform

    ! Constrained (omega, alpha_g, beta_g, alp_nig) -> unconstrained p(4).
    subroutine garch_nig_inv_transform(omega, alpha_g, beta_g, alp, p)
        real(dp), intent(in)  :: omega, alpha_g, beta_g, alp
        real(dp), intent(out) :: p(4)
        call garch_inv_transform(omega, alpha_g, beta_g, p(1:3))
        p(4) = log((alp - alp_min) / (alp_max - alp))
    end subroutine garch_nig_inv_transform

    ! Simulate GARCH(1,1) with symmetric NIG(alp) innovations.
    ! epsilon = sqrt(V) * Z,  V ~ InvGauss(mu=1, lambda=alp^2),  Z ~ N(0,1).
    subroutine garch_nig_simulate(omega, alpha_g, beta_g, alp, n, seed_val, y)
        real(dp), intent(in)  :: omega, alpha_g, beta_g, alp
        integer,  intent(in)  :: n, seed_val
        real(dp), intent(out) :: y(n)
        real(dp) :: h, eps
        integer  :: i, sz
        integer, allocatable :: seed_arr(:)

        call random_seed(size=sz)
        allocate(seed_arr(sz))
        seed_arr = seed_val
        call random_seed(put=seed_arr)
        deallocate(seed_arr)

        h = omega / (1.0_dp - alpha_g - beta_g)
        do i = 1, n
            eps = random_nig_sym(alp)
            y(i) = sqrt(h) * eps
            h    = omega + alpha_g * y(i)**2 + beta_g * h
        end do
    end subroutine garch_nig_simulate

    ! Negative log-likelihood (normalised by n) and gradient w.r.t. unconstrained p(4).
    subroutine garch_nig_obj(p, np, f, g)
        integer,  intent(in)  :: np
        real(dp), intent(in)  :: p(np)
        real(dp), intent(out) :: f
        real(dp), intent(out) :: g(np)

        real(dp) :: omega, alpha_g, beta_g, alp
        real(dp) :: h, dh_dom, dh_dal, dh_dbe
        real(dp) :: grad_om, grad_al, grad_be, grad_ap
        real(dp) :: h_unc, u, u2, a2pu2, r, lk1, ratio, factor, dalpd4
        integer  :: t

        call garch_nig_transform(p, omega, alpha_g, beta_g, alp)
        ! d(alp)/d(p4) = (alp - alp_min) * (alp_max - alp) / alp_rng
        dalpd4 = (alp - alp_min) * (alp_max - alp) / alp_rng

        h_unc  = omega / (1.0_dp - alpha_g - beta_g)
        h      = h_unc
        dh_dom = 1.0_dp / (1.0_dp - alpha_g - beta_g)
        dh_dal = h_unc  / (1.0_dp - alpha_g - beta_g)
        dh_dbe = h_unc  / (1.0_dp - alpha_g - beta_g)

        ! constant part: n*(log(pi) - 2*log(alp) - alp^2)
        f       = real(gni_nobs, dp) * (log(pi) - 2.0_dp*log(alp) - alp*alp)
        grad_om = 0.0_dp
        grad_al = 0.0_dp
        grad_be = 0.0_dp
        ! constant contribution to d(NLL)/d(alp) from n*(-2*log(alp) - alp^2):
        ! d/d(alp)[n*(-2*log(alp) - alp^2)] = n*(-2/alp - 2*alp);
        ! the +n/alp from Σ Q/r is folded in below (see nig_bessel derivation), giving n*(-1/alp - 2*alp).
        grad_ap = real(gni_nobs, dp) * (-1.0_dp/alp - 2.0_dp*alp)

        do t = 1, gni_nobs
            u      = gni_obs(t) / sqrt(h)
            u2     = u * u
            a2pu2  = alp*alp + u2           ! alp^2 + u^2
            r      = alp * sqrt(a2pu2)      ! alp * sqrt(alp^2 + u^2)
            call bessel_k01(r, lk1, ratio)  ! log K1(r) and K0(r)/K1(r)

            f       = f + 0.5_dp*log(h) + 0.5_dp*log(a2pu2) - lk1

            ! gradient w.r.t. h
            factor  = (1.0_dp - u2*(2.0_dp + r*ratio)/a2pu2) / (2.0_dp*h)
            grad_om = grad_om + factor * dh_dom
            grad_al = grad_al + factor * dh_dal
            grad_be = grad_be + factor * dh_dbe

            ! per-observation contribution to d(NLL)/d(alp):
            ! from d/d(alp)[0.5*log(alp^2+u^2) - log K1(r)].
            ! The 1/r * Q term from d/d(alp)[-log K1] = 1/alp + alp/(alp^2+u^2),
            ! the 1/alp part combines with the n*(-2/alp) constant to give n*(-1/alp-2*alp) above.
            ! What remains per obs: K0/K1*Q + alp/(alp^2+u^2) + alp/(alp^2+u^2)  = K0/K1*Q + 2*alp/(alp^2+u^2)
            grad_ap = grad_ap + ratio*(2.0_dp*alp*alp + u2)/sqrt(a2pu2) + 2.0_dp*alp/a2pu2

            ! advance recurrences to h_{t+1}
            dh_dom = 1.0_dp             + beta_g * dh_dom
            dh_dal = gni_obs(t)**2      + beta_g * dh_dal
            dh_dbe = h                  + beta_g * dh_dbe
            h      = omega + alpha_g * gni_obs(t)**2 + beta_g * h
        end do

        ! chain rule: constrained gradients -> unconstrained
        g(1) =  grad_om * omega
        g(2) =  grad_al * alpha_g*(1.0_dp - alpha_g) - grad_be * alpha_g*beta_g
        g(3) = -grad_al * alpha_g*beta_g             + grad_be * beta_g*(1.0_dp - beta_g)
        g(4) =  grad_ap * dalpd4

        f = f / gni_nobs
        g = g / gni_nobs
    end subroutine garch_nig_obj

end module garch_nig_mod

! GARCH(1,1) module: simulation, parameter transforms, negative log-likelihood
! with analytical gradients with respect to unconstrained parameters.
!
! Parameterisation (constrained):
!   h_t = omega + alpha*y_{t-1}^2 + beta*h_{t-1},  h_1 = omega/(1-alpha-beta)
!   omega > 0,  alpha > 0,  beta > 0,  alpha+beta < 1
!
! Unconstrained map  p -> (omega, alpha, beta):
!   omega = exp(p1)
!   alpha = exp(p2) / (1 + exp(p2) + exp(p3))
!   beta  = exp(p3) / (1 + exp(p2) + exp(p3))
! Inverse:
!   p1 = log(omega)
!   p2 = log(alpha / (1-alpha-beta))
!   p3 = log(beta  / (1-alpha-beta))

module garch_mod
    use kind_mod,       only: dp
    use math_const_mod, only: log_sqrt_2pi, two_pi
    implicit none
    private

    ! data shared with objective callback
    real(dp), allocatable, save :: g_obs(:)
    integer,               save :: g_nobs = 0
    real(dp),              save :: g_target_var = 0.0_dp   ! variance target for VT variant

    public :: garch_set_data, garch_simulate, garch_obj, &
              garch_transform, garch_inv_transform, &
              garch_vt_set_target, garch_vt_transform, garch_vt_inv_transform, garch_vt_obj

contains

    ! Store observations for use by garch_obj.
    subroutine garch_set_data(y, n)
        integer,  intent(in) :: n
        real(dp), intent(in) :: y(n)
        if (allocated(g_obs)) deallocate(g_obs)
        allocate(g_obs(n))
        g_obs  = y
        g_nobs = n
    end subroutine garch_set_data

    ! Unconstrained p(3) -> constrained (omega, alpha, beta).
    subroutine garch_transform(p, omega, alpha, beta)
        real(dp), intent(in)  :: p(3)
        real(dp), intent(out) :: omega, alpha, beta
        real(dp) :: e2, e3, s
        omega = exp(p(1))
        e2    = exp(p(2))
        e3    = exp(p(3))
        s     = 1.0_dp + e2 + e3
        alpha = e2 / s
        beta  = e3 / s
    end subroutine garch_transform

    ! Constrained (omega, alpha, beta) -> unconstrained p(3).
    ! Requires omega > 0, alpha > 0, beta > 0, alpha+beta < 1.
    subroutine garch_inv_transform(omega, alpha, beta, p)
        real(dp), intent(in)  :: omega, alpha, beta
        real(dp), intent(out) :: p(3)
        real(dp) :: gamma
        gamma = 1.0_dp - alpha - beta
        p(1)  = log(omega)
        p(2)  = log(alpha / gamma)
        p(3)  = log(beta  / gamma)
    end subroutine garch_inv_transform

    ! Simulate n observations from GARCH(1,1) using Box-Muller normals.
    ! h_1 is set to the unconditional variance omega/(1-alpha-beta).
    subroutine garch_simulate(omega, alpha, beta, n, seed_val, y)
        real(dp), intent(in)  :: omega, alpha, beta
        integer,  intent(in)  :: n, seed_val
        real(dp), intent(out) :: y(n)
        real(dp) :: h, u1, u2, eps
        integer  :: i, sz
        integer, allocatable :: seed_arr(:)

        call random_seed(size=sz)
        allocate(seed_arr(sz))
        seed_arr = seed_val
        call random_seed(put=seed_arr)
        deallocate(seed_arr)

        h = omega / (1.0_dp - alpha - beta)
        do i = 1, n
            ! guard against log(0) in Box-Muller
            do
                call random_number(u1)
                if (u1 > 0.0_dp) exit
            end do
            call random_number(u2)
            eps  = sqrt(-2.0_dp * log(u1)) * cos(two_pi * u2)
            y(i) = sqrt(h) * eps
            h    = omega + alpha * y(i)**2 + beta * h
        end do
    end subroutine garch_simulate

    ! Negative log-likelihood and its gradient w.r.t. unconstrained p.
    !
    ! f = -loglik = 0.5 * sum_t [ log(h_t) + y_t^2/h_t ]   (constant omitted)
    !
    ! Recurrences for dh_t/d_theta (theta = omega, alpha, beta):
    !   dh_t/d_omega = 1            + beta * dh_{t-1}/d_omega
    !   dh_t/d_alpha = y_{t-1}^2   + beta * dh_{t-1}/d_alpha
    !   dh_t/d_beta  = h_{t-1}     + beta * dh_{t-1}/d_beta
    !   initialised from h_1 = omega/(1-alpha-beta)
    !
    ! Chain rule from constrained to unconstrained (softmax Jacobian):
    !   d_omega/d_p1 = omega
    !   d_alpha/d_p2 = alpha*(1-alpha),  d_alpha/d_p3 = -alpha*beta
    !   d_beta /d_p2 = -alpha*beta,      d_beta /d_p3 =  beta*(1-beta)
    subroutine garch_obj(p, np, f, g)
        integer,  intent(in)  :: np
        real(dp), intent(in)  :: p(np)
        real(dp), intent(out) :: f
        real(dp), intent(out) :: g(np)

        real(dp) :: omega, alpha, beta
        real(dp) :: h, dh_dom, dh_dal, dh_dbe
        real(dp) :: grad_om, grad_al, grad_be
        real(dp) :: h_unc, factor
        integer  :: t

        call garch_transform(p, omega, alpha, beta)

        h_unc  = omega / (1.0_dp - alpha - beta)
        h      = h_unc
        dh_dom = 1.0_dp / (1.0_dp - alpha - beta)
        dh_dal = h_unc  / (1.0_dp - alpha - beta)
        dh_dbe = h_unc  / (1.0_dp - alpha - beta)

        f       = real(g_nobs, dp) * log_sqrt_2pi
        grad_om = 0.0_dp
        grad_al = 0.0_dp
        grad_be = 0.0_dp

        do t = 1, g_nobs
            factor  =  1.0_dp/h - g_obs(t)**2 / h**2
            f       = f       + 0.5_dp * (log(h) + g_obs(t)**2 / h)
            grad_om = grad_om + 0.5_dp * factor * dh_dom
            grad_al = grad_al + 0.5_dp * factor * dh_dal
            grad_be = grad_be + 0.5_dp * factor * dh_dbe
            ! advance recurrences to h_{t+1}
            dh_dom = 1.0_dp          + beta * dh_dom
            dh_dal = g_obs(t)**2     + beta * dh_dal
            dh_dbe = h               + beta * dh_dbe
            h      = omega + alpha * g_obs(t)**2 + beta * h
        end do

        ! chain rule: constrained gradients -> unconstrained gradients
        g(1) =  grad_om * omega
        g(2) =  grad_al * alpha*(1.0_dp - alpha) - grad_be * alpha*beta
        g(3) = -grad_al * alpha*beta             + grad_be * beta*(1.0_dp - beta)

        ! normalise by n so that f and ||g|| are O(1) regardless of sample size
        f = f / g_nobs
        g = g / g_nobs
    end subroutine garch_obj

    ! ── Variance-targeted GARCH ───────────────────────────────────────────────
    ! omega is derived: omega = target*(1-alpha-beta), so h_1 = target exactly.
    ! Reduces stage-1 parameters from 3 to 2 per asset.

    subroutine garch_vt_set_target(target_var)
        ! Store per-asset unconditional variance target before calling garch_vt_obj.
        real(dp), intent(in) :: target_var  ! sample variance of demeaned returns
        g_target_var = target_var
    end subroutine garch_vt_set_target

    subroutine garch_vt_transform(p, alpha, beta)
        ! p(1:2) -> (alpha, beta) via 2-component softmax.
        real(dp), intent(in)  :: p(2)
        real(dp), intent(out) :: alpha, beta
        real(dp) :: e1, e2, s
        e1    = exp(p(1));  e2 = exp(p(2))
        s     = 1.0_dp + e1 + e2
        alpha = e1 / s;     beta = e2 / s
    end subroutine garch_vt_transform

    subroutine garch_vt_inv_transform(alpha, beta, p)
        ! (alpha, beta) -> p(1:2).
        real(dp), intent(in)  :: alpha, beta
        real(dp), intent(out) :: p(2)
        real(dp) :: gamma
        gamma = 1.0_dp - alpha - beta
        p(1)  = log(alpha / gamma)
        p(2)  = log(beta  / gamma)
    end subroutine garch_vt_inv_transform

    subroutine garch_vt_obj(p, np, f, g)
        ! NLL/n and gradient for variance-targeted GARCH(1,1). np=2.
        !
        ! omega = target*(1-alpha-beta),  h_1 = target  (dh_1/d* = 0).
        !
        ! Derivative recurrences (domega/dalpha = -target, domega/dbeta = -target):
        !   dh_t/dalpha = -target + y_{t-1}^2 + beta * dh_{t-1}/dalpha
        !   dh_t/dbeta  = -target + h_{t-1}   + beta * dh_{t-1}/dbeta
        !
        ! Chain rule (2-component softmax Jacobian):
        !   g(1) =  grad_al * alpha*(1-alpha) - grad_be * alpha*beta
        !   g(2) = -grad_al * alpha*beta      + grad_be * beta*(1-beta)
        integer,  intent(in)  :: np
        real(dp), intent(in)  :: p(np)
        real(dp), intent(out) :: f, g(np)
        real(dp) :: alpha, beta, omega, h, dh_dal, dh_dbe, grad_al, grad_be, factor
        integer  :: t
        call garch_vt_transform(p, alpha, beta)
        omega  = g_target_var * (1.0_dp - alpha - beta)
        h      = g_target_var    ! h_1 = target by construction
        dh_dal = 0.0_dp          ! dh_1/dalpha = 0
        dh_dbe = 0.0_dp          ! dh_1/dbeta  = 0
        f       = real(g_nobs, dp) * log_sqrt_2pi
        grad_al = 0.0_dp;  grad_be = 0.0_dp
        do t = 1, g_nobs
            factor  = 1.0_dp/h - g_obs(t)**2 / h**2
            f       = f       + 0.5_dp * (log(h) + g_obs(t)**2 / h)
            grad_al = grad_al + 0.5_dp * factor * dh_dal
            grad_be = grad_be + 0.5_dp * factor * dh_dbe
            dh_dal  = -g_target_var + g_obs(t)**2 + beta * dh_dal
            dh_dbe  = -g_target_var + h            + beta * dh_dbe
            h       = omega + alpha * g_obs(t)**2  + beta * h
        end do
        g(1) =  grad_al * alpha*(1.0_dp - alpha) - grad_be * alpha*beta
        g(2) = -grad_al * alpha*beta             + grad_be * beta*(1.0_dp - beta)
        f = f / g_nobs;  g = g / g_nobs
    end subroutine garch_vt_obj

end module garch_mod

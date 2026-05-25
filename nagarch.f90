! GARCH(1,1) with NAGARCH variance equation and Normal innovations.
!
! Variance equation:
!   h_t = omega + alpha*q_{t-1}^2 + beta*h_{t-1}
!   q_t = y_t - theta*sqrt(h_t) by default, or max(theta*sqrt(h_t)-y_t, 0)
!         when zero-above-shift news impact is enabled.
!   h_1 = omega / (1 - alpha*E(q_t^2/h_t) - beta)
!   omega > 0,  alpha > 0,  beta > 0,  theta unconstrained
!   Stationarity: alpha*(1+theta^2) + beta < 1
!
! Unconstrained parameterisation  p(4) -> (omega, alpha, beta, theta):
!   Let aa = alpha*(1+theta^2),  S = 1 + exp(p2) + exp(p3)
!   omega = exp(p1)
!   aa    = exp(p2) / S       ! satisfies aa > 0, aa+beta < 1
!   beta  = exp(p3) / S
!   alpha = aa / (1+theta^2)
!   theta = p4                ! unconstrained
! Inverse:
!   p1 = log(omega)
!   p2 = log(alpha*(1+theta^2) / (1 - alpha*(1+theta^2) - beta))
!   p3 = log(beta              / (1 - alpha*(1+theta^2) - beta))
!   p4 = theta

module nagarch_mod
    use kind_mod,       only: dp
    use math_const_mod, only: log_sqrt_2pi, two_pi, sqrt2
    implicit none
    private

    real(dp), allocatable, save :: gna_obs(:)  ! stored observations
    integer,               save :: gna_nobs = 0
    real(dp),              save :: gna_target_var = 0.0_dp   ! variance target for VT variant
    logical,               save :: gna_zero_above_shift = .false.

    public :: nagarch_set_data, nagarch_simulate, nagarch_obj, &
              nagarch_transform, nagarch_inv_transform, &
              nagarch_vt_set_target, nagarch_vt_transform, nagarch_vt_inv_transform, nagarch_vt_obj, &
              nagarch_set_news_impact, nagarch_news_impact_zero_above

contains

    subroutine nagarch_set_news_impact(zero_above_shift)
        logical, intent(in) :: zero_above_shift
        gna_zero_above_shift = zero_above_shift
    end subroutine nagarch_set_news_impact

    logical function nagarch_news_impact_zero_above()
        nagarch_news_impact_zero_above = gna_zero_above_shift
    end function nagarch_news_impact_zero_above

    subroutine nagarch_set_data(y, n)
        ! Store observations for use by nagarch_obj.
        integer,  intent(in) :: n    ! number of observations
        real(dp), intent(in) :: y(n) ! return series
        if (allocated(gna_obs)) deallocate(gna_obs)
        allocate(gna_obs(n))
        gna_obs  = y
        gna_nobs = n
    end subroutine nagarch_set_data

    subroutine nagarch_transform(p, omega, alpha, beta, theta)
        ! Unconstrained p(4) -> constrained (omega, alpha, beta, theta).
        real(dp), intent(in)  :: p(4)
        real(dp), intent(out) :: omega, alpha, beta, theta  ! NAGARCH parameters
        real(dp) :: e2, e3, s, aa, shift_moment, dmoment_unused
        omega = exp(p(1))
        theta = p(4)
        call nagarch_shift_moments(theta, shift_moment, dmoment_unused)
        e2    = exp(p(2))
        e3    = exp(p(3))
        s     = 1.0_dp + e2 + e3
        aa    = e2 / s
        beta  = e3 / s
        alpha = aa / shift_moment
    end subroutine nagarch_transform

    subroutine nagarch_inv_transform(omega, alpha, beta, theta, p)
        ! Constrained (omega, alpha, beta, theta) -> unconstrained p(4).
        ! Requires omega > 0, alpha > 0, beta > 0, alpha*(1+theta^2)+beta < 1.
        real(dp), intent(in)  :: omega, alpha, beta, theta  ! NAGARCH parameters
        real(dp), intent(out) :: p(4)
        real(dp) :: aa, gam, shift_moment, dmoment_unused
        call nagarch_shift_moments(theta, shift_moment, dmoment_unused)
        aa   = alpha * shift_moment
        gam  = 1.0_dp - aa - beta            ! stationarity slack
        p(1) = log(omega)
        p(2) = log(aa   / gam)
        p(3) = log(beta / gam)
        p(4) = theta
    end subroutine nagarch_inv_transform

    subroutine nagarch_simulate(omega, alpha, beta, theta, n, seed_val, y)
        ! Simulate n observations from NAGARCH(1,1) with Normal innovations.
        ! Uses Box-Muller; h_1 is set to the unconditional variance.
        real(dp), intent(in)  :: omega, alpha, beta, theta  ! NAGARCH parameters
        integer,  intent(in)  :: n                          ! number of observations
        integer,  intent(in)  :: seed_val                   ! RNG seed
        real(dp), intent(out) :: y(n)                       ! simulated return series
        real(dp) :: h, sqrth, r, u1, u2, eps, shift_moment, dmoment_unused
        integer  :: i, sz
        integer, allocatable :: seed_arr(:)

        call random_seed(size=sz)
        allocate(seed_arr(sz))
        seed_arr = seed_val
        call random_seed(put=seed_arr)
        deallocate(seed_arr)

        call nagarch_shift_moments(theta, shift_moment, dmoment_unused)
        h = omega / (1.0_dp - alpha*shift_moment - beta)
        do i = 1, n
            do
                call random_number(u1)
                if (u1 > 0.0_dp) exit
            end do
            call random_number(u2)
            eps   = sqrt(-2.0_dp * log(u1)) * cos(two_pi * u2)
            sqrth = sqrt(h)
            y(i)  = sqrth * eps
            r     = y(i) - theta * sqrth     ! shifted residual
            if (gna_zero_above_shift .and. r > 0.0_dp) r = 0.0_dp
            h     = omega + alpha*r**2 + beta*h
        end do
    end subroutine nagarch_simulate

    subroutine nagarch_obj(p, np, f, g)
        ! NLL and gradient w.r.t. unconstrained p for NAGARCH(1,1) with Normal noise.
        !
        ! f = (1/n) * [ n*log_sqrt_2pi + 0.5*sum_t(log(h_t) + y_t^2/h_t) ]
        !
        ! Define r_t = y_{t-1} - theta*sqrt(h_{t-1})
        !        kappa_t = beta - alpha*theta*r_t / sqrt(h_{t-1})
        !
        ! Variance derivative recurrences (forcing terms differ per parameter):
        !   dh_t/d_omega = 1                          + kappa_t * dh_{t-1}/d_omega
        !   dh_t/d_alpha = r_t^2                      + kappa_t * dh_{t-1}/d_alpha
        !   dh_t/d_beta  = h_{t-1}                    + kappa_t * dh_{t-1}/d_beta
        !   dh_t/d_theta = -2*alpha*r_t*sqrt(h_{t-1}) + kappa_t * dh_{t-1}/d_theta
        !
        ! Initial values from h_1 = omega/D, D = 1 - alpha*(1+theta^2) - beta:
        !   dh_1/d_omega = 1/D
        !   dh_1/d_alpha = h_1*(1+theta^2)/D
        !   dh_1/d_beta  = h_1/D
        !   dh_1/d_theta = 2*alpha*theta*h_1/D
        !
        ! Chain rule to unconstrained params (aa = alpha*(1+theta^2), s2 = 1+theta^2):
        !   g(1) =  grad_om * omega
        !   g(2) =  grad_al * alpha*(1-aa) - grad_be * aa*beta
        !   g(3) = -grad_al * alpha*beta   + grad_be * beta*(1-beta)
        !   g(4) =  grad_al * (-2*theta*alpha/s2) + grad_th
        integer,  intent(in)  :: np     ! number of parameters (must be 4)
        real(dp), intent(in)  :: p(np)  ! unconstrained parameters
        real(dp), intent(out) :: f      ! NLL/n
        real(dp), intent(out) :: g(np)  ! gradient of NLL/n

        real(dp) :: omega, alpha, beta, theta
        real(dp) :: h, sqrth, r, kappa, s2, D, h_unc, aa, dmom
        real(dp) :: dh_dom, dh_dal, dh_dbe, dh_dth
        real(dp) :: grad_om, grad_al, grad_be, grad_th, factor
        logical  :: active
        integer  :: t

        call nagarch_transform(p, omega, alpha, beta, theta)

        call nagarch_shift_moments(theta, s2, dmom)
        D     = 1.0_dp - alpha*s2 - beta
        h_unc = omega / D

        ! initial variance derivatives from h_1 = omega/D
        h      = h_unc
        dh_dom = 1.0_dp / D
        dh_dal = h_unc * s2 / D
        dh_dbe = h_unc / D
        dh_dth = alpha * dmom * h_unc / D

        f       = real(gna_nobs, dp) * log_sqrt_2pi
        grad_om = 0.0_dp
        grad_al = 0.0_dp
        grad_be = 0.0_dp
        grad_th = 0.0_dp

        do t = 1, gna_nobs
            factor  =  1.0_dp/h - gna_obs(t)**2 / h**2
            f       = f       + 0.5_dp * (log(h) + gna_obs(t)**2 / h)
            grad_om = grad_om + 0.5_dp * factor * dh_dom
            grad_al = grad_al + 0.5_dp * factor * dh_dal
            grad_be = grad_be + 0.5_dp * factor * dh_dbe
            grad_th = grad_th + 0.5_dp * factor * dh_dth
            ! update recurrences for h_{t+1} using r_{t+1} = y_t - theta*sqrt(h_t)
            sqrth  = sqrt(h)
            r      = gna_obs(t) - theta * sqrth
            active = .true.
            if (gna_zero_above_shift .and. r > 0.0_dp) then
                r = 0.0_dp
                active = .false.
            end if
            kappa  = beta
            if (active) kappa = beta - alpha * theta * r / sqrth
            dh_dom = 1.0_dp                + kappa * dh_dom
            dh_dal = r**2                  + kappa * dh_dal
            dh_dbe = h                     + kappa * dh_dbe
            dh_dth = kappa * dh_dth
            if (active) dh_dth = dh_dth - 2.0_dp*alpha*r*sqrth
            h      = omega + alpha*r**2 + beta*h
        end do

        ! chain rule: constrained gradients -> unconstrained gradients
        aa   = alpha * s2
        g(1) =  grad_om * omega
        g(2) =  grad_al * alpha*(1.0_dp - aa) - grad_be * aa*beta
        g(3) = -grad_al * alpha*beta           + grad_be * beta*(1.0_dp - beta)
        g(4) =  grad_al * (-alpha*dmom/s2) + grad_th

        f = f / gna_nobs
        g = g / gna_nobs
    end subroutine nagarch_obj

    ! ── Variance-targeted NAGARCH ─────────────────────────────────────────────
    ! omega = target*(1-alpha*(1+theta^2)-beta), h_1 = target (dh_1/d* = 0).
    ! Reduces stage-1 parameters from 4 to 3 per asset.
    !
    ! Free parameters: p(1) = logit(aa), p(2) = logit(beta), p(3) = theta
    ! where aa = alpha*(1+theta^2), same softmax as standard but without omega.
    !
    ! Derivative recurrences (kappa = beta - alpha*theta*r_t/sqrt(h_t)):
    !   dh_{t+1}/dalpha = -target*s2 + r_t^2              + kappa*dh_t/dalpha
    !   dh_{t+1}/dbeta  = -target    + h_t                + kappa*dh_t/dbeta
    !   dh_{t+1}/dtheta = -2*target*alpha*theta - 2*alpha*r_t*sqrt(h_t) + kappa*dh_t/dtheta
    ! where s2 = 1+theta^2.

    subroutine nagarch_vt_set_target(target_var)
        ! Store per-asset unconditional variance target before calling nagarch_vt_obj.
        real(dp), intent(in) :: target_var  ! sample variance of demeaned returns
        gna_target_var = target_var
    end subroutine nagarch_vt_set_target

    subroutine nagarch_vt_transform(p, alpha, beta, theta)
        ! p(1:3) -> (alpha, beta, theta); aa=alpha*(1+theta^2) via softmax on p(1:2).
        real(dp), intent(in)  :: p(3)
        real(dp), intent(out) :: alpha, beta, theta
        real(dp) :: e1, e2, s, aa, shift_moment, dmoment_unused
        theta = p(3)
        call nagarch_shift_moments(theta, shift_moment, dmoment_unused)
        e1    = exp(p(1));  e2 = exp(p(2))
        s     = 1.0_dp + e1 + e2
        aa    = e1 / s;     beta  = e2 / s
        alpha = aa / shift_moment
    end subroutine nagarch_vt_transform

    subroutine nagarch_vt_inv_transform(alpha, beta, theta, p)
        ! (alpha, beta, theta) -> p(1:3).
        real(dp), intent(in)  :: alpha, beta, theta
        real(dp), intent(out) :: p(3)
        real(dp) :: aa, gam, shift_moment, dmoment_unused
        call nagarch_shift_moments(theta, shift_moment, dmoment_unused)
        aa   = alpha * shift_moment
        gam  = 1.0_dp - aa - beta
        p(1) = log(aa   / gam)
        p(2) = log(beta / gam)
        p(3) = theta
    end subroutine nagarch_vt_inv_transform

    subroutine nagarch_vt_obj(p, np, f, g)
        ! NLL/n and gradient for variance-targeted NAGARCH(1,1). np=3.
        !
        ! Chain rule to unconstrained (aa = alpha*(1+theta^2)):
        !   g(1) =  grad_al * alpha*(1-aa)       - grad_be * aa*beta
        !   g(2) = -grad_al * alpha*beta          + grad_be * beta*(1-beta)
        !   g(3) =  grad_al * (-2*theta*alpha/s2) + grad_th
        integer,  intent(in)  :: np
        real(dp), intent(in)  :: p(np)
        real(dp), intent(out) :: f, g(np)
        real(dp) :: alpha, beta, theta, omega, s2, aa, dmom
        real(dp) :: h, sqrth, r, kappa
        real(dp) :: dh_dal, dh_dbe, dh_dth
        real(dp) :: grad_al, grad_be, grad_th, factor
        logical  :: active
        integer  :: t
        call nagarch_vt_transform(p, alpha, beta, theta)
        call nagarch_shift_moments(theta, s2, dmom)
        omega = gna_target_var * (1.0_dp - alpha*s2 - beta)
        h     = gna_target_var     ! h_1 = target by construction
        dh_dal = 0.0_dp;  dh_dbe = 0.0_dp;  dh_dth = 0.0_dp
        f       = real(gna_nobs, dp) * log_sqrt_2pi
        grad_al = 0.0_dp;  grad_be = 0.0_dp;  grad_th = 0.0_dp
        do t = 1, gna_nobs
            factor  = 1.0_dp/h - gna_obs(t)**2 / h**2
            f       = f       + 0.5_dp * (log(h) + gna_obs(t)**2 / h)
            grad_al = grad_al + 0.5_dp * factor * dh_dal
            grad_be = grad_be + 0.5_dp * factor * dh_dbe
            grad_th = grad_th + 0.5_dp * factor * dh_dth
            sqrth  = sqrt(h)
            r      = gna_obs(t) - theta * sqrth
            active = .true.
            if (gna_zero_above_shift .and. r > 0.0_dp) then
                r = 0.0_dp
                active = .false.
            end if
            kappa  = beta
            if (active) kappa = beta - alpha * theta * r / sqrth
            dh_dal = -gna_target_var*s2                     + r**2             + kappa * dh_dal
            dh_dbe = -gna_target_var                        + h                + kappa * dh_dbe
            dh_dth = -gna_target_var*alpha*dmom + kappa * dh_dth
            if (active) dh_dth = dh_dth - 2.0_dp*alpha*r*sqrth
            h      = omega + alpha*r**2 + beta*h
        end do
        aa   = alpha * s2
        g(1) =  grad_al * alpha*(1.0_dp - aa)        - grad_be * aa*beta
        g(2) = -grad_al * alpha*beta                 + grad_be * beta*(1.0_dp - beta)
        g(3) =  grad_al * (-alpha*dmom/s2)            + grad_th
        f = f / gna_nobs;  g = g / gna_nobs
    end subroutine nagarch_vt_obj

    subroutine nagarch_shift_moments(theta, moment, dmoment)
        real(dp), intent(in)  :: theta
        real(dp), intent(out) :: moment, dmoment
        real(dp) :: cdf, pdf
        if (gna_zero_above_shift) then
            cdf     = 0.5_dp * (1.0_dp + erf(theta / sqrt2))
            pdf     = exp(-0.5_dp * theta**2 - log_sqrt_2pi)
            moment  = (1.0_dp + theta**2) * cdf + theta * pdf
            dmoment = 2.0_dp * (theta * cdf + pdf)
        else
            moment  = 1.0_dp + theta**2
            dmoment = 2.0_dp * theta
        end if
    end subroutine nagarch_shift_moments

end module nagarch_mod

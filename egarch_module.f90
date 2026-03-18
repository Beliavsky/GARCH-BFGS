! EGARCH(1,1) with Normal innovations (Nelson, 1991).
!
! Variance equation (log form):
!   log(h_t) = omega + beta*log(h_{t-1}) + alpha*(|z_{t-1}| - c) + gamma*z_{t-1}
!   where z_t = y_t/sqrt(h_t),  c = sqrt(2/pi)  (E[|z|] under Normal)
!   log(h_1) = omega / (1 - beta)
!   omega, alpha, gamma unconstrained;  |beta| < 1
!
! Unconstrained parameterisation  p(4) -> (omega, alpha, gamma, beta):
!   omega = p(1)                        (free)
!   alpha = p(2)                        (free)
!   gamma = p(3)                        (free)
!   beta  = tanh(p(4))                  (enforces |beta| < 1)
! Inverse:
!   p(1) = omega
!   p(2) = alpha
!   p(3) = gamma
!   p(4) = arctanh(beta) = 0.5*log((1+beta)/(1-beta))

module egarch_module
    use kind_mod,       only: dp
    use math_const_mod, only: log_sqrt_2pi, two_pi, pi
    implicit none
    private

    real(dp), allocatable, save :: ega_obs(:)  ! stored observations
    integer,               save :: ega_nobs = 0

    public :: egarch_set_data, egarch_simulate, egarch_obj, &
              egarch_transform, egarch_inv_transform

contains

    subroutine egarch_set_data(y, n)
        ! Store observations for use by egarch_obj.
        integer,  intent(in) :: n    ! number of observations
        real(dp), intent(in) :: y(n) ! return series
        if (allocated(ega_obs)) deallocate(ega_obs)
        allocate(ega_obs(n))
        ega_obs  = y
        ega_nobs = n
    end subroutine egarch_set_data

    subroutine egarch_transform(p, omega, alpha, gamma, beta)
        ! Unconstrained p(4) -> constrained (omega, alpha, gamma, beta).
        real(dp), intent(in)  :: p(4)
        real(dp), intent(out) :: omega, alpha, gamma, beta
        omega = p(1)
        alpha = p(2)
        gamma = p(3)
        beta  = tanh(p(4))
    end subroutine egarch_transform

    subroutine egarch_inv_transform(omega, alpha, gamma, beta, p)
        ! Constrained (omega, alpha, gamma, beta) -> unconstrained p(4).
        ! Requires |beta| < 1.
        real(dp), intent(in)  :: omega, alpha, gamma, beta
        real(dp), intent(out) :: p(4)
        p(1) = omega
        p(2) = alpha
        p(3) = gamma
        p(4) = 0.5_dp * log((1.0_dp + beta) / (1.0_dp - beta))  ! arctanh
    end subroutine egarch_inv_transform

    subroutine egarch_simulate(omega, alpha, gamma, beta, n, seed_val, y)
        ! Simulate n observations from EGARCH(1,1) with Normal innovations.
        ! Uses Box-Muller; log(h_1) is set to the unconditional log-variance.
        real(dp), intent(in)  :: omega, alpha, gamma, beta
        integer,  intent(in)  :: n, seed_val
        real(dp), intent(out) :: y(n)
        real(dp) :: lh, h, z, u1, u2, eps, c_eg
        integer  :: i, sz
        integer, allocatable :: seed_arr(:)

        call random_seed(size=sz)
        allocate(seed_arr(sz))
        seed_arr = seed_val
        call random_seed(put=seed_arr)
        deallocate(seed_arr)

        c_eg = sqrt(2.0_dp / pi)
        lh   = omega / (1.0_dp - beta)   ! unconditional log-variance
        do i = 1, n
            h = exp(lh)
            do
                call random_number(u1)
                if (u1 > 0.0_dp) exit
            end do
            call random_number(u2)
            eps  = sqrt(-2.0_dp * log(u1)) * cos(two_pi * u2)
            y(i) = sqrt(h) * eps
            z    = eps                    ! z_t = y_t/sqrt(h_t) = eps
            lh   = omega + beta*lh + alpha*(abs(z) - c_eg) + gamma*z
        end do
    end subroutine egarch_simulate

    subroutine egarch_obj(p, np, f, g)
        ! NLL and gradient w.r.t. unconstrained p for EGARCH(1,1) with Normal noise.
        !
        ! f = (1/n) * [ n*log_sqrt_2pi + 0.5*sum_t(lh_t + y_t^2*exp(-lh_t)) ]
        !
        ! Define z_t = y_t/sqrt(h_t),  c = sqrt(2/pi)
        !        kappa_t = beta - 0.5*(alpha*|z_{t-1}| + gamma*z_{t-1})
        !          (from d(lh_t)/d(lh_{t-1}) via dz_{t-1}/d(lh_{t-1}) = -z_{t-1}/2)
        !
        ! Log-variance derivative recurrences:
        !   dlh_t/d_omega = 1                    + kappa_t * dlh_{t-1}/d_omega
        !   dlh_t/d_alpha = (|z_{t-1}| - c)      + kappa_t * dlh_{t-1}/d_alpha
        !   dlh_t/d_gamma = z_{t-1}              + kappa_t * dlh_{t-1}/d_gamma
        !   dlh_t/d_beta  = lh_{t-1}             + kappa_t * dlh_{t-1}/d_beta
        !
        ! Initial values from lh_1 = omega/(1-beta):
        !   dlh_1/d_omega = 1/(1-beta)
        !   dlh_1/d_alpha = 0
        !   dlh_1/d_gamma = 0
        !   dlh_1/d_beta  = omega/(1-beta)^2 = lh_1/(1-beta)
        !
        ! Chain rule to unconstrained params:
        !   g(1) = grad_om                       (omega = p1, free)
        !   g(2) = grad_al                       (alpha = p2, free)
        !   g(3) = grad_ga                       (gamma = p3, free)
        !   g(4) = grad_be * (1 - beta^2)        (beta = tanh(p4))
        integer,  intent(in)  :: np
        real(dp), intent(in)  :: p(np)
        real(dp), intent(out) :: f
        real(dp), intent(out) :: g(np)

        real(dp) :: omega, alpha, gamma, beta
        real(dp) :: lh, h, z, abs_z, kappa, c_eg
        real(dp) :: dlh_dom, dlh_dal, dlh_dga, dlh_dbe
        real(dp) :: grad_om, grad_al, grad_ga, grad_be, factor
        integer  :: t

        call egarch_transform(p, omega, alpha, gamma, beta)

        c_eg = sqrt(2.0_dp / pi)
        lh   = omega / (1.0_dp - beta)   ! unconditional log-variance

        ! initial log-variance derivatives from lh_1 = omega/(1-beta)
        dlh_dom = 1.0_dp / (1.0_dp - beta)
        dlh_dal = 0.0_dp
        dlh_dga = 0.0_dp
        dlh_dbe = lh / (1.0_dp - beta)   ! = omega/(1-beta)^2

        f       = real(ega_nobs, dp) * log_sqrt_2pi
        grad_om = 0.0_dp
        grad_al = 0.0_dp
        grad_ga = 0.0_dp
        grad_be = 0.0_dp

        do t = 1, ega_nobs
            h       = exp(lh)
            z       = ega_obs(t) / sqrt(h)
            factor  = 0.5_dp * (1.0_dp - z**2)
            f       = f       + 0.5_dp * (lh + ega_obs(t)**2 / h)
            grad_om = grad_om + factor * dlh_dom
            grad_al = grad_al + factor * dlh_dal
            grad_ga = grad_ga + factor * dlh_dga
            grad_be = grad_be + factor * dlh_dbe
            ! update recurrences for lh_{t+1}
            abs_z   = abs(z)
            kappa   = beta - 0.5_dp * (alpha * abs_z + gamma * z)
            dlh_dom = 1.0_dp             + kappa * dlh_dom
            dlh_dal = (abs_z - c_eg)     + kappa * dlh_dal
            dlh_dga = z                  + kappa * dlh_dga
            dlh_dbe = lh                 + kappa * dlh_dbe
            lh      = omega + beta*lh + alpha*(abs_z - c_eg) + gamma*z
        end do

        ! chain rule: constrained gradients -> unconstrained gradients
        g(1) = grad_om
        g(2) = grad_al
        g(3) = grad_ga
        g(4) = grad_be * (1.0_dp - beta**2)

        f = f / ega_nobs
        g = g / ega_nobs
    end subroutine egarch_obj

end module egarch_module

! GARCH(1,1) with GJR-GARCH variance equation and Normal innovations.
!
! Variance equation:
!   h_t = omega + alpha*y_{t-1}^2 + gamma*I(y_{t-1}<0)*y_{t-1}^2 + beta*h_{t-1}
!   h_1 = omega / (1 - alpha - gamma/2 - beta)
!   omega > 0,  alpha > 0,  gamma >= 0,  beta > 0
!   Stationarity: alpha + gamma/2 + beta < 1
!   (E[I(z<0)*z^2] = 1/2 under Normal, so effective persistence = alpha+gamma/2+beta)
!
! Unconstrained parameterisation  p(4) -> (omega, alpha, gamma, beta):
!   Let S = 1 + exp(p2) + exp(p3) + exp(p4)
!   omega = exp(p1)
!   alpha = exp(p2) / S
!   gamma = 2*exp(p3) / S     (gamma/2 enters the softmax, ensuring gamma >= 0)
!   beta  = exp(p4) / S
!   Stationarity: alpha + gamma/2 + beta = (exp(p2)+exp(p3)+exp(p4))/S < 1
! Inverse:
!   p1 = log(omega)
!   p2 = log(alpha        / slack)    where slack = 1 - alpha - gamma/2 - beta
!   p3 = log(gamma/(2*slack))
!   p4 = log(beta         / slack)

module gjr_module
    use kind_mod,       only: dp
    use math_const_mod, only: log_sqrt_2pi, two_pi
    implicit none
    private

    real(dp), allocatable, save :: gjr_obs(:)  ! stored observations
    integer,               save :: gjr_nobs = 0

    public :: gjr_set_data, gjr_simulate, gjr_obj, &
              gjr_transform, gjr_inv_transform

contains

    subroutine gjr_set_data(y, n)
        ! Store observations for use by gjr_obj.
        integer,  intent(in) :: n    ! number of observations
        real(dp), intent(in) :: y(n) ! return series
        if (allocated(gjr_obs)) deallocate(gjr_obs)
        allocate(gjr_obs(n))
        gjr_obs  = y
        gjr_nobs = n
    end subroutine gjr_set_data

    subroutine gjr_transform(p, omega, alpha, gamma, beta)
        ! Unconstrained p(4) -> constrained (omega, alpha, gamma, beta).
        real(dp), intent(in)  :: p(4)
        real(dp), intent(out) :: omega, alpha, gamma, beta  ! GJR parameters
        real(dp) :: e2, e3, e4, s
        omega = exp(p(1))
        e2    = exp(p(2))
        e3    = exp(p(3))
        e4    = exp(p(4))
        s     = 1.0_dp + e2 + e3 + e4
        alpha = e2 / s
        gamma = 2.0_dp * e3 / s   ! gamma/2 in the softmax
        beta  = e4 / s
    end subroutine gjr_transform

    subroutine gjr_inv_transform(omega, alpha, gamma, beta, p)
        ! Constrained (omega, alpha, gamma, beta) -> unconstrained p(4).
        ! Requires omega > 0, alpha > 0, gamma >= 0, beta > 0,
        !          alpha + gamma/2 + beta < 1.
        real(dp), intent(in)  :: omega, alpha, gamma, beta  ! GJR parameters
        real(dp), intent(out) :: p(4)
        real(dp) :: slack
        slack = 1.0_dp - alpha - 0.5_dp*gamma - beta   ! stationarity slack
        p(1)  = log(omega)
        p(2)  = log(alpha / slack)
        p(3)  = log(0.5_dp*gamma / slack)
        p(4)  = log(beta  / slack)
    end subroutine gjr_inv_transform

    subroutine gjr_simulate(omega, alpha, gamma, beta, n, seed_val, y)
        ! Simulate n observations from GJR-GARCH(1,1) with Normal innovations.
        ! Uses Box-Muller; h_1 is set to the unconditional variance.
        real(dp), intent(in)  :: omega, alpha, gamma, beta  ! GJR parameters
        integer,  intent(in)  :: n                          ! number of observations
        integer,  intent(in)  :: seed_val                   ! RNG seed
        real(dp), intent(out) :: y(n)                       ! simulated return series
        real(dp) :: h, ind, u1, u2, eps
        integer  :: i, sz
        integer, allocatable :: seed_arr(:)

        call random_seed(size=sz)
        allocate(seed_arr(sz))
        seed_arr = seed_val
        call random_seed(put=seed_arr)
        deallocate(seed_arr)

        h = omega / (1.0_dp - alpha - 0.5_dp*gamma - beta)
        do i = 1, n
            do
                call random_number(u1)
                if (u1 > 0.0_dp) exit
            end do
            call random_number(u2)
            eps  = sqrt(-2.0_dp * log(u1)) * cos(two_pi * u2)
            y(i) = sqrt(h) * eps
            if (y(i) < 0.0_dp) then
                ind = 1.0_dp
            else
                ind = 0.0_dp
            end if
            h = omega + (alpha + gamma*ind)*y(i)**2 + beta*h
        end do
    end subroutine gjr_simulate

    subroutine gjr_obj(p, np, f, g)
        ! NLL and gradient w.r.t. unconstrained p for GJR-GARCH(1,1) with Normal noise.
        !
        ! f = (1/n) * [ n*log_sqrt_2pi + 0.5*sum_t(log(h_t) + y_t^2/h_t) ]
        !
        ! Variance derivative recurrences (multiplier is beta throughout):
        !   dh_t/d_omega = 1                  + beta * dh_{t-1}/d_omega
        !   dh_t/d_alpha = y_{t-1}^2          + beta * dh_{t-1}/d_alpha
        !   dh_t/d_gamma = I_{t-1}*y_{t-1}^2 + beta * dh_{t-1}/d_gamma
        !   dh_t/d_beta  = h_{t-1}            + beta * dh_{t-1}/d_beta
        !
        ! Initial values from h_1 = omega/D, D = 1 - alpha - gamma/2 - beta:
        !   dh_1/d_omega = 1/D
        !   dh_1/d_alpha = h_1/D
        !   dh_1/d_gamma = h_1/(2*D)     (factor 1/2: d(D)/d(gamma) = -1/2)
        !   dh_1/d_beta  = h_1/D
        !
        ! Chain rule to unconstrained params (g2 = gamma/2):
        !   g(1) =  grad_om * omega
        !   g(2) =  alpha * (grad_al*(1-alpha) - 2*grad_ga*g2 - grad_be*beta)
        !   g(3) =  g2    * (-grad_al*alpha + 2*grad_ga*(1-g2) - grad_be*beta)
        !   g(4) =  beta  * (-grad_al*alpha - 2*grad_ga*g2 + grad_be*(1-beta))
        integer,  intent(in)  :: np     ! number of parameters (must be 4)
        real(dp), intent(in)  :: p(np)  ! unconstrained parameters
        real(dp), intent(out) :: f      ! NLL/n
        real(dp), intent(out) :: g(np)  ! gradient of NLL/n

        real(dp) :: omega, alpha, gamma, beta
        real(dp) :: h, ind, D, h_unc, g2
        real(dp) :: dh_dom, dh_dal, dh_dga, dh_dbe
        real(dp) :: grad_om, grad_al, grad_ga, grad_be, factor
        integer  :: t

        call gjr_transform(p, omega, alpha, gamma, beta)

        g2    = 0.5_dp * gamma
        D     = 1.0_dp - alpha - g2 - beta
        h_unc = omega / D

        ! initial variance derivatives from h_1 = omega/D
        h      = h_unc
        dh_dom = 1.0_dp / D
        dh_dal = h_unc / D
        dh_dga = h_unc / (2.0_dp * D)   ! factor 1/2 from d(D)/d(gamma) = -1/2
        dh_dbe = h_unc / D

        f       = real(gjr_nobs, dp) * log_sqrt_2pi
        grad_om = 0.0_dp
        grad_al = 0.0_dp
        grad_ga = 0.0_dp
        grad_be = 0.0_dp

        do t = 1, gjr_nobs
            factor  =  1.0_dp/h - gjr_obs(t)**2 / h**2
            f       = f       + 0.5_dp * (log(h) + gjr_obs(t)**2 / h)
            grad_om = grad_om + 0.5_dp * factor * dh_dom
            grad_al = grad_al + 0.5_dp * factor * dh_dal
            grad_ga = grad_ga + 0.5_dp * factor * dh_dga
            grad_be = grad_be + 0.5_dp * factor * dh_dbe
            ! update recurrences for h_{t+1}
            if (gjr_obs(t) < 0.0_dp) then
                ind = 1.0_dp
            else
                ind = 0.0_dp
            end if
            dh_dom = 1.0_dp              + beta * dh_dom
            dh_dal = gjr_obs(t)**2       + beta * dh_dal
            dh_dga = ind*gjr_obs(t)**2   + beta * dh_dga
            dh_dbe = h                   + beta * dh_dbe
            h      = omega + (alpha + gamma*ind)*gjr_obs(t)**2 + beta*h
        end do

        ! chain rule: constrained gradients -> unconstrained gradients
        g(1) =  grad_om * omega
        g(2) =  alpha * (grad_al*(1.0_dp-alpha) - 2.0_dp*grad_ga*g2 - grad_be*beta)
        g(3) =  g2    * (-grad_al*alpha + 2.0_dp*grad_ga*(1.0_dp-g2) - grad_be*beta)
        g(4) =  beta  * (-grad_al*alpha - 2.0_dp*grad_ga*g2 + grad_be*(1.0_dp-beta))

        f = f / gjr_nobs
        g = g / gjr_nobs
    end subroutine gjr_obj

end module gjr_module

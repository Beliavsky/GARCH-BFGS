! GARCH(1,1) with standardised logistic innovations.
!
! Innovation model:  y_t = sqrt(h_t) * epsilon_t
!   epsilon_t ~ Logistic(0, sqrt(3)/pi): PDF f(x) = c*exp(-c*x)/(1+exp(-c*x))^2
!   c = pi/sqrt(3)  ensures unit variance.
!
! Negative log-likelihood per observation:
!   nll_t = 0.5*log(h_t) - log(c) + |u_t| + 2*log(1 + exp(-|u_t|))
!   u_t = c * y_t / sqrt(h_t)
!
! Numerically stable form uses log(1+exp(-|u|)) which is always O(1) since |u|>=0.
!
! Gradient w.r.t. h_t:
!   d(nll_t)/dh_t = (1 - |u_t| * tanh(|u_t|/2)) / (2*h_t)
!
! Simulation via inverse CDF of logistic:
!   epsilon = (1/c) * log(U/(1-U)),  U ~ Uniform(0,1)

module garch_logistic_module
    use kind_mod,       only: dp
    use math_const_mod, only: pi, sqrt3
    use garch_module,   only: garch_transform, garch_inv_transform
    implicit none
    private

    real(dp), parameter :: c_logis = pi / sqrt3  ! scale c = pi/sqrt(3) for unit variance

    real(dp), allocatable, save :: glo_obs(:)
    integer,               save :: glo_nobs = 0

    public :: garch_logistic_set_data, garch_logistic_simulate, garch_logistic_obj

contains

    subroutine garch_logistic_set_data(y, n)
        integer,  intent(in) :: n
        real(dp), intent(in) :: y(n)
        if (allocated(glo_obs)) deallocate(glo_obs)
        allocate(glo_obs(n))
        glo_obs  = y
        glo_nobs = n
    end subroutine garch_logistic_set_data

    ! Simulate GARCH(1,1) with standardised logistic innovations via inverse CDF.
    ! epsilon = (1/c) * log(U/(1-U)),  U ~ Uniform(0,1)
    subroutine garch_logistic_simulate(omega, alpha, beta, n, seed_val, y)
        real(dp), intent(in)  :: omega, alpha, beta
        integer,  intent(in)  :: n, seed_val
        real(dp), intent(out) :: y(n)
        real(dp) :: h, u, eps
        integer  :: i, sz
        integer, allocatable :: seed_arr(:)

        call random_seed(size=sz)
        allocate(seed_arr(sz))
        seed_arr = seed_val
        call random_seed(put=seed_arr)
        deallocate(seed_arr)

        h = omega / (1.0_dp - alpha - beta)
        do i = 1, n
            do
                call random_number(u)
                if (u > 0.0_dp .and. u < 1.0_dp) exit
            end do
            eps  = log(u / (1.0_dp - u)) / c_logis
            y(i) = sqrt(h) * eps
            h    = omega + alpha * y(i)**2 + beta * h
        end do
    end subroutine garch_logistic_simulate

    ! Negative log-likelihood (normalised by n) and gradient w.r.t. unconstrained p(3).
    subroutine garch_logistic_obj(p, np, f, g)
        integer,  intent(in)  :: np
        real(dp), intent(in)  :: p(np)
        real(dp), intent(out) :: f
        real(dp), intent(out) :: g(np)

        real(dp) :: omega, alpha, beta
        real(dp) :: h, dh_dom, dh_dal, dh_dbe
        real(dp) :: grad_om, grad_al, grad_be
        real(dp) :: h_unc, u, au, factor
        integer  :: t

        call garch_transform(p, omega, alpha, beta)

        h_unc  = omega / (1.0_dp - alpha - beta)
        h      = h_unc
        dh_dom = 1.0_dp / (1.0_dp - alpha - beta)
        dh_dal = h_unc  / (1.0_dp - alpha - beta)
        dh_dbe = h_unc  / (1.0_dp - alpha - beta)

        f       = -real(glo_nobs, dp) * log(c_logis)
        grad_om = 0.0_dp
        grad_al = 0.0_dp
        grad_be = 0.0_dp

        do t = 1, glo_nobs
            u      = c_logis * glo_obs(t) / sqrt(h)
            au     = abs(u)
            f      = f + 0.5_dp*log(h) + au + 2.0_dp*log(1.0_dp + exp(-au))
            factor = (1.0_dp - au*tanh(0.5_dp*au)) / (2.0_dp*h)
            grad_om = grad_om + factor * dh_dom
            grad_al = grad_al + factor * dh_dal
            grad_be = grad_be + factor * dh_dbe
            ! advance recurrences to h_{t+1}
            dh_dom = 1.0_dp           + beta * dh_dom
            dh_dal = glo_obs(t)**2    + beta * dh_dal
            dh_dbe = h                + beta * dh_dbe
            h      = omega + alpha * glo_obs(t)**2 + beta * h
        end do

        ! chain rule: constrained gradients -> unconstrained (softmax Jacobian)
        g(1) =  grad_om * omega
        g(2) =  grad_al * alpha*(1.0_dp - alpha) - grad_be * alpha*beta
        g(3) = -grad_al * alpha*beta             + grad_be * beta*(1.0_dp - beta)

        f = f / glo_nobs
        g = g / glo_nobs
    end subroutine garch_logistic_obj

end module garch_logistic_module

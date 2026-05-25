! GARCH(1,1) with standardised hyperbolic-secant innovations.
!
! Innovation model:  y_t = sqrt(h_t) * epsilon_t
!   epsilon_t ~ sech: PDF f(x) = (1/2)*sech(pi*x/2), mean 0, variance 1.
!
! Conditional negative log-likelihood per observation:
!   nll_t = 0.5*log(h_t) + log(2) + log(cosh(u_t)),  u_t = pi*y_t / (2*sqrt(h_t))
!
! Gradient w.r.t. h_t:
!   d(nll_t)/d(h_t) = (1 - u_t*tanh(u_t)) / (2*h_t)
!
! Simulation via inverse CDF of sech distribution:
!   epsilon = (2/pi) * log(tan(pi*U/2)),  U ~ Uniform(0,1)
!
! Unconstrained parameterisation is identical to Normal GARCH: p(3) -> (omega, alpha, beta).

module garch_sech_mod
    use kind_mod,       only: dp
    use math_const_mod, only: pi, log2
    use garch_mod,   only: garch_transform, garch_inv_transform
    implicit none
    private

    real(dp), allocatable, save :: gs_obs(:)
    integer,               save :: gs_nobs = 0

    public :: garch_sech_set_data, garch_sech_simulate, garch_sech_obj

contains

    subroutine garch_sech_set_data(y, n)
        integer,  intent(in) :: n
        real(dp), intent(in) :: y(n)
        if (allocated(gs_obs)) deallocate(gs_obs)
        allocate(gs_obs(n))
        gs_obs  = y
        gs_nobs = n
    end subroutine garch_sech_set_data

    ! Simulate GARCH(1,1) with standardised sech innovations via inverse CDF.
    subroutine garch_sech_simulate(omega, alpha, beta, n, seed_val, y)
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
            eps  = (2.0_dp / pi) * log(tan(pi * u / 2.0_dp))
            y(i) = sqrt(h) * eps
            h    = omega + alpha * y(i)**2 + beta * h
        end do
    end subroutine garch_sech_simulate

    ! Negative log-likelihood (normalised by n) and gradient w.r.t. unconstrained p(3).
    subroutine garch_sech_obj(p, np, f, g)
        integer,  intent(in)  :: np
        real(dp), intent(in)  :: p(np)
        real(dp), intent(out) :: f
        real(dp), intent(out) :: g(np)

        real(dp) :: omega, alpha, beta
        real(dp) :: h, dh_dom, dh_dal, dh_dbe
        real(dp) :: grad_om, grad_al, grad_be
        real(dp) :: h_unc, u, factor
        integer  :: t

        call garch_transform(p, omega, alpha, beta)

        h_unc  = omega / (1.0_dp - alpha - beta)
        h      = h_unc
        dh_dom = 1.0_dp / (1.0_dp - alpha - beta)
        dh_dal = h_unc  / (1.0_dp - alpha - beta)
        dh_dbe = h_unc  / (1.0_dp - alpha - beta)

        f       = real(gs_nobs, dp) * log2
        grad_om = 0.0_dp
        grad_al = 0.0_dp
        grad_be = 0.0_dp

        do t = 1, gs_nobs
            u       = pi * gs_obs(t) / (2.0_dp * sqrt(h))
            f       = f + 0.5_dp*log(h) + log_cosh(u)
            factor  = (1.0_dp - u*tanh(u)) / (2.0_dp * h)
            grad_om = grad_om + factor * dh_dom
            grad_al = grad_al + factor * dh_dal
            grad_be = grad_be + factor * dh_dbe
            ! advance recurrences to h_{t+1}
            dh_dom = 1.0_dp          + beta * dh_dom
            dh_dal = gs_obs(t)**2    + beta * dh_dal
            dh_dbe = h               + beta * dh_dbe
            h      = omega + alpha * gs_obs(t)**2 + beta * h
        end do

        ! chain rule: constrained gradients -> unconstrained (softmax Jacobian)
        g(1) =  grad_om * omega
        g(2) =  grad_al * alpha*(1.0_dp - alpha) - grad_be * alpha*beta
        g(3) = -grad_al * alpha*beta             + grad_be * beta*(1.0_dp - beta)

        f = f / gs_nobs
        g = g / gs_nobs
    end subroutine garch_sech_obj

    ! Numerically stable log(cosh(u)): avoids overflow for large |u|.
    ! log(cosh(u)) = |u| + log(1 + exp(-2|u|)) - log(2)
    pure function log_cosh(u) result(lc)
        real(dp), intent(in) :: u
        real(dp) :: lc, au
        au = abs(u)
        lc = au + log(1.0_dp + exp(-2.0_dp*au)) - log2
    end function log_cosh

end module garch_sech_mod

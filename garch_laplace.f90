! GARCH(1,1) with standardised Laplace (double-exponential) innovations.
!
! Innovation model:  y_t = sqrt(h_t) * epsilon_t
!   epsilon_t ~ Laplace(0, 1/sqrt(2)): PDF f(x) = (1/sqrt(2))*exp(-sqrt(2)*|x|)
!   This gives mean 0 and variance 1.
!   Equivalent to GED with nu=1 but with simpler arithmetic (no Gamma calls).
!
! Negative log-likelihood per observation:
!   nll_t = 0.5*log(h_t) + sqrt(2)*|y_t|/sqrt(h_t) + 0.5*log(2)
!
! Gradient w.r.t. h_t:
!   d(nll_t)/dh_t = (1 - sqrt(2)*|y_t|/sqrt(h_t)) / (2*h_t)
!
! Simulation via difference of unit exponentials:
!   epsilon = (E1 - E2) / sqrt(2),  E1,E2 ~ Exp(1) iid  =>  variance = 1

module garch_laplace_mod
    use kind_mod,       only: dp
    use math_const_mod, only: sqrt2, half_log2
    use garch_mod,   only: garch_transform, garch_inv_transform
    implicit none
    private

    real(dp), allocatable, save :: gl_obs(:)
    integer,               save :: gl_nobs = 0

    public :: garch_laplace_set_data, garch_laplace_simulate, garch_laplace_obj

contains

    subroutine garch_laplace_set_data(y, n)
        integer,  intent(in) :: n
        real(dp), intent(in) :: y(n)
        if (allocated(gl_obs)) deallocate(gl_obs)
        allocate(gl_obs(n))
        gl_obs  = y
        gl_nobs = n
    end subroutine garch_laplace_set_data

    ! Simulate GARCH(1,1) with standardised Laplace innovations.
    ! epsilon = (E1 - E2) / sqrt(2),  E1,E2 ~ Exp(1) via -log(U).
    subroutine garch_laplace_simulate(omega, alpha, beta, n, seed_val, y)
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
            do
                call random_number(u1)
                if (u1 > 0.0_dp) exit
            end do
            do
                call random_number(u2)
                if (u2 > 0.0_dp) exit
            end do
            eps  = (-log(u1) + log(u2)) / sqrt2
            y(i) = sqrt(h) * eps
            h    = omega + alpha * y(i)**2 + beta * h
        end do
    end subroutine garch_laplace_simulate

    ! Negative log-likelihood (normalised by n) and gradient w.r.t. unconstrained p(3).
    subroutine garch_laplace_obj(p, np, f, g)
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

        f       = real(gl_nobs, dp) * half_log2
        grad_om = 0.0_dp
        grad_al = 0.0_dp
        grad_be = 0.0_dp

        do t = 1, gl_nobs
            f       = f + 0.5_dp*log(h) + sqrt2*abs(gl_obs(t))/sqrt(h)
            factor  = (1.0_dp - sqrt2*abs(gl_obs(t))/sqrt(h)) / (2.0_dp*h)
            grad_om = grad_om + factor * dh_dom
            grad_al = grad_al + factor * dh_dal
            grad_be = grad_be + factor * dh_dbe
            ! advance recurrences to h_{t+1}
            dh_dom = 1.0_dp           + beta * dh_dom
            dh_dal = gl_obs(t)**2     + beta * dh_dal
            dh_dbe = h                + beta * dh_dbe
            h      = omega + alpha * gl_obs(t)**2 + beta * h
        end do

        ! chain rule: constrained gradients -> unconstrained (softmax Jacobian)
        g(1) =  grad_om * omega
        g(2) =  grad_al * alpha*(1.0_dp - alpha) - grad_be * alpha*beta
        g(3) = -grad_al * alpha*beta             + grad_be * beta*(1.0_dp - beta)

        f = f / gl_nobs
        g = g / gl_nobs
    end subroutine garch_laplace_obj

end module garch_laplace_mod

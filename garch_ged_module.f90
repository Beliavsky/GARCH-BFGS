! GARCH(1,1) with Generalised Error Distribution (GED) innovations.
!
! Innovation model:  y_t = sqrt(h_t) * epsilon_t
!   epsilon_t ~ standardised GED(nu): PDF f(x) = nu/(2^(1+1/nu)*lam*Gamma(1/nu)) * exp(-1/2*|x/lam|^nu)
!   lam(nu) = 2^(-1/nu) * sqrt(Gamma(1/nu)/Gamma(3/nu))  ensures unit variance
!   nu=2 -> Normal,  nu=1 -> Laplace,  nu<2 heavy tails,  nu>2 light tails
!
! Negative log-likelihood per observation:
!   nll_t = 0.5*log(h_t) + 0.5*a_t^nu + f_const(nu)
!   a_t   = |y_t| / (lam(nu) * sqrt(h_t))
!   f_const = -log(nu) + log(2) + 1.5*log_gamma(1/nu) - 0.5*log_gamma(3/nu)
!
! Gradient w.r.t. h_t:
!   d(nll_t)/dh_t = (1 - (nu/2)*a_t^nu) / (2*h_t)
!
! Gradient w.r.t. nu (per observation):
!   d(nll_t)/dnu = dfc/dnu + 0.5*a_t^nu*(log(a_t) - c_nu)
!   dfc/dnu = -1/nu + (3/(2*nu^2)) * [digamma(3/nu) - digamma(1/nu)]
!   c_nu    = [log(2) - 0.5*digamma(1/nu) + 1.5*digamma(3/nu)] / nu
!
! Unconstrained parameterisation:
!   p(1:3) -> (omega, alpha, beta)  via softmax (same as Normal GARCH)
!   p(4)   -> nu = exp(p(4))        ensures nu > 0;  d(nu)/d(p4) = nu

module garch_ged_module
    use kind_mod,          only: dp
    use garch_module,      only: garch_transform, garch_inv_transform
    use random_mod,        only: random_gamma
    use distributions_mod, only: ged_lambda
    use special_mod,       only: digamma
    implicit none
    private

    real(dp), allocatable, save :: gg_obs(:)
    integer,               save :: gg_nobs = 0

    public :: garch_ged_set_data, garch_ged_simulate, garch_ged_obj, &
              garch_ged_transform, garch_ged_inv_transform

contains

    subroutine garch_ged_set_data(y, n)
        integer,  intent(in) :: n
        real(dp), intent(in) :: y(n)
        if (allocated(gg_obs)) deallocate(gg_obs)
        allocate(gg_obs(n))
        gg_obs  = y
        gg_nobs = n
    end subroutine garch_ged_set_data

    ! Unconstrained p(4) -> constrained (omega, alpha, beta, nu).
    subroutine garch_ged_transform(p, omega, alpha, beta, nu)
        real(dp), intent(in)  :: p(4)
        real(dp), intent(out) :: omega, alpha, beta, nu
        call garch_transform(p(1:3), omega, alpha, beta)
        nu = exp(p(4))
    end subroutine garch_ged_transform

    ! Constrained (omega, alpha, beta, nu) -> unconstrained p(4).
    subroutine garch_ged_inv_transform(omega, alpha, beta, nu, p)
        real(dp), intent(in)  :: omega, alpha, beta, nu
        real(dp), intent(out) :: p(4)
        call garch_inv_transform(omega, alpha, beta, p(1:3))
        p(4) = log(nu)
    end subroutine garch_ged_inv_transform

    ! Simulate GARCH(1,1) with standardised GED(nu) innovations.
    ! Uses:  epsilon = sign * lam(nu) * (2*G)^(1/nu),  G ~ Gamma(1/nu, 1)
    subroutine garch_ged_simulate(omega, alpha, beta, nu, n, seed_val, y)
        real(dp), intent(in)  :: omega, alpha, beta, nu
        integer,  intent(in)  :: n, seed_val
        real(dp), intent(out) :: y(n)
        real(dp) :: h, lam, g, u, eps
        integer  :: i, sz
        integer, allocatable :: seed_arr(:)

        call random_seed(size=sz)
        allocate(seed_arr(sz))
        seed_arr = seed_val
        call random_seed(put=seed_arr)
        deallocate(seed_arr)

        lam = ged_lambda(nu)
        h   = omega / (1.0_dp - alpha - beta)
        do i = 1, n
            g = random_gamma(1.0_dp/nu)
            call random_number(u)
            eps  = lam * (2.0_dp*g)**(1.0_dp/nu) * sign(1.0_dp, u - 0.5_dp)
            y(i) = sqrt(h) * eps
            h    = omega + alpha * y(i)**2 + beta * h
        end do
    end subroutine garch_ged_simulate

    ! Negative log-likelihood (normalised by n) and gradient w.r.t. unconstrained p(4).
    subroutine garch_ged_obj(p, np, f, g)
        integer,  intent(in)  :: np
        real(dp), intent(in)  :: p(np)
        real(dp), intent(out) :: f
        real(dp), intent(out) :: g(np)

        real(dp) :: omega, alpha, beta, nu
        real(dp) :: h, dh_dom, dh_dal, dh_dbe
        real(dp) :: grad_om, grad_al, grad_be, grad_nu
        real(dp) :: h_unc, lam, a, a_nu, factor, f_const, dfc_dnu, c_nu
        integer  :: t

        call garch_ged_transform(p, omega, alpha, beta, nu)

        lam     = ged_lambda(nu)
        f_const = -log(nu) + log(2.0_dp) &
                  + 1.5_dp*log_gamma(1.0_dp/nu) - 0.5_dp*log_gamma(3.0_dp/nu)
        dfc_dnu = -1.0_dp/nu &
                  + (1.5_dp/nu**2) * (digamma(3.0_dp/nu) - digamma(1.0_dp/nu))
        c_nu    = (log(2.0_dp) - 0.5_dp*digamma(1.0_dp/nu) &
                  + 1.5_dp*digamma(3.0_dp/nu)) / nu

        h_unc  = omega / (1.0_dp - alpha - beta)
        h      = h_unc
        dh_dom = 1.0_dp / (1.0_dp - alpha - beta)
        dh_dal = h_unc  / (1.0_dp - alpha - beta)
        dh_dbe = h_unc  / (1.0_dp - alpha - beta)

        f       = real(gg_nobs, dp) * f_const
        grad_om = 0.0_dp
        grad_al = 0.0_dp
        grad_be = 0.0_dp
        grad_nu = 0.0_dp

        do t = 1, gg_nobs
            a      = abs(gg_obs(t)) / (lam * sqrt(h))
            a_nu   = a**nu
            f      = f + 0.5_dp*log(h) + 0.5_dp*a_nu
            factor = (1.0_dp - 0.5_dp*nu*a_nu) / (2.0_dp*h)
            grad_om = grad_om + factor * dh_dom
            grad_al = grad_al + factor * dh_dal
            grad_be = grad_be + factor * dh_dbe
            if (a > 0.0_dp) then
                grad_nu = grad_nu + dfc_dnu + 0.5_dp*a_nu*(log(a) - c_nu)
            else
                grad_nu = grad_nu + dfc_dnu
            end if
            ! advance recurrences to h_{t+1}
            dh_dom = 1.0_dp           + beta * dh_dom
            dh_dal = gg_obs(t)**2     + beta * dh_dal
            dh_dbe = h                + beta * dh_dbe
            h      = omega + alpha * gg_obs(t)**2 + beta * h
        end do

        ! chain rule: constrained gradients -> unconstrained (softmax + exp Jacobians)
        g(1) =  grad_om * omega
        g(2) =  grad_al * alpha*(1.0_dp - alpha) - grad_be * alpha*beta
        g(3) = -grad_al * alpha*beta             + grad_be * beta*(1.0_dp - beta)
        g(4) =  grad_nu * nu

        f = f / gg_nobs
        g = g / gg_nobs
    end subroutine garch_ged_obj

    ! ---- private routines ----------------------------------------

end module garch_ged_module

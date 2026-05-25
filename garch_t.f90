! GARCH(1,1) with standardised Student-t innovations.
!
! Innovation model:  y_t = sqrt(h_t) * epsilon_t
!   epsilon_t ~ standardised t(nu): zero mean, unit variance (requires nu > 2)
!   epsilon_t = z_t * sqrt((nu-2)/nu),  z_t ~ standard t(nu)
!
! Conditional negative log-likelihood per observation:
!   f_t = -log_gamma((nu+1)/2) + log_gamma(nu/2) + 0.5*log(pi*(nu-2))
!         + 0.5*log(h_t) + 0.5*(nu+1)*log(1 + y_t^2/((nu-2)*h_t))
!
! Unconstrained parameterisation  p(4) -> (omega, alpha, beta, nu):
!   p(1:3) -> (omega, alpha, beta)  via softmax map (same as garch_mod)
!   p(4)   -> nu = 2 + 98/(1+exp(-p(4)))   ensures 2 < nu < 100
! Inverse:
!   p(4) = log((nu-2)/(100-nu))
! Jacobian:
!   dnu/dp4 = (nu-2)*(100-nu)/98

module garch_t_mod
    use kind_mod,       only: dp
    use math_const_mod, only: pi
    use garch_mod,   only: garch_transform, garch_inv_transform
    use random_mod,     only: random_normal, random_gamma, random_t_std
    use special_mod,    only: digamma
    implicit none
    private

    real(dp), allocatable, save :: gt_obs(:)
    integer,               save :: gt_nobs = 0

    public :: garch_t_set_data, garch_t_simulate, garch_t_obj, &
              garch_t_transform, garch_t_inv_transform

contains

    subroutine garch_t_set_data(y, n)
        integer,  intent(in) :: n
        real(dp), intent(in) :: y(n)
        if (allocated(gt_obs)) deallocate(gt_obs)
        allocate(gt_obs(n))
        gt_obs  = y
        gt_nobs = n
    end subroutine garch_t_set_data

    ! Unconstrained p(4) -> constrained (omega, alpha, beta, nu).
    subroutine garch_t_transform(p, omega, alpha, beta, nu)
        real(dp), intent(in)  :: p(4)
        real(dp), intent(out) :: omega, alpha, beta, nu
        call garch_transform(p(1:3), omega, alpha, beta)
        nu = 2.0_dp + 98.0_dp / (1.0_dp + exp(-p(4)))
    end subroutine garch_t_transform

    ! Constrained (omega, alpha, beta, nu) -> unconstrained p(4).
    subroutine garch_t_inv_transform(omega, alpha, beta, nu, p)
        real(dp), intent(in)  :: omega, alpha, beta, nu
        real(dp), intent(out) :: p(4)
        call garch_inv_transform(omega, alpha, beta, p(1:3))
        p(4) = log((nu - 2.0_dp) / (100.0_dp - nu))
    end subroutine garch_t_inv_transform

    ! Simulate GARCH(1,1) with standardised t(nu) innovations.
    ! h_1 is set to the unconditional variance omega/(1-alpha-beta).
    subroutine garch_t_simulate(omega, alpha, beta, nu, n, seed_val, y)
        real(dp), intent(in)  :: omega, alpha, beta, nu
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

        h = omega / (1.0_dp - alpha - beta)
        do i = 1, n
            eps = random_t_std(nu)
            y(i) = sqrt(h) * eps
            h    = omega + alpha * y(i)**2 + beta * h
        end do
    end subroutine garch_t_simulate

    ! Negative log-likelihood (normalised by n) and gradient w.r.t. p(4).
    !
    ! Gradient w.r.t. h_t:
    !   df_t/dh_t = 1/(2h_t) * (1 - (nu+1)*u_t/(1+u_t)),  u_t = y_t^2/((nu-2)*h_t)
    !
    ! Gradient w.r.t. nu (per observation):
    !   df_t/dnu = -0.5*psi((nu+1)/2) + 0.5*psi(nu/2) + 0.5/(nu-2)
    !              + 0.5*log(1+u_t) - 0.5*(nu+1)*u_t/((nu-2)*(1+u_t))
    !   where psi = digamma
    !
    ! Chain rule to unconstrained:
    !   g(1:3) same softmax Jacobian as normal GARCH
    !   g(4) = sum_t df_t/dnu * (nu-2)   [d_nu/d_p4 = exp(p4) = nu-2]
    subroutine garch_t_obj(p, np, f, g)
        integer,  intent(in)  :: np
        real(dp), intent(in)  :: p(np)
        real(dp), intent(out) :: f
        real(dp), intent(out) :: g(np)

        real(dp) :: omega, alpha, beta, nu
        real(dp) :: h, dh_dom, dh_dal, dh_dbe
        real(dp) :: grad_om, grad_al, grad_be, grad_nu
        real(dp) :: h_unc, u, log1pu, factor, f_const, dnu_const
        integer  :: t

        call garch_t_transform(p, omega, alpha, beta, nu)

        ! constant part of f_t and df_t/dnu (independent of t)
        f_const   = -log_gamma(0.5_dp*(nu+1.0_dp)) + log_gamma(0.5_dp*nu) &
                    + 0.5_dp*log(pi*(nu-2.0_dp))
        dnu_const = -0.5_dp*digamma(0.5_dp*(nu+1.0_dp)) &
                    + 0.5_dp*digamma(0.5_dp*nu)          &
                    + 0.5_dp/(nu-2.0_dp)

        h_unc  = omega / (1.0_dp - alpha - beta)
        h      = h_unc
        dh_dom = 1.0_dp / (1.0_dp - alpha - beta)
        dh_dal = h_unc  / (1.0_dp - alpha - beta)
        dh_dbe = h_unc  / (1.0_dp - alpha - beta)

        f       = real(gt_nobs, dp) * f_const
        grad_om = 0.0_dp
        grad_al = 0.0_dp
        grad_be = 0.0_dp
        grad_nu = 0.0_dp

        do t = 1, gt_nobs
            u       = gt_obs(t)**2 / ((nu-2.0_dp)*h)
            log1pu  = log(1.0_dp + u)

            f       = f       + 0.5_dp*log(h) + 0.5_dp*(nu+1.0_dp)*log1pu
            factor  = 0.5_dp/h * (1.0_dp - (nu+1.0_dp)*u/(1.0_dp+u))
            grad_om = grad_om + factor * dh_dom
            grad_al = grad_al + factor * dh_dal
            grad_be = grad_be + factor * dh_dbe
            grad_nu = grad_nu + dnu_const &
                              + 0.5_dp*log1pu &
                              - 0.5_dp*(nu+1.0_dp)*u/((nu-2.0_dp)*(1.0_dp+u))

            ! advance recurrences to h_{t+1}
            dh_dom = 1.0_dp          + beta * dh_dom
            dh_dal = gt_obs(t)**2    + beta * dh_dal
            dh_dbe = h               + beta * dh_dbe
            h      = omega + alpha * gt_obs(t)**2 + beta * h
        end do

        ! chain rule: constrained gradients -> unconstrained
        g(1) =  grad_om * omega
        g(2) =  grad_al * alpha*(1.0_dp - alpha) - grad_be * alpha*beta
        g(3) = -grad_al * alpha*beta             + grad_be * beta*(1.0_dp - beta)
        g(4) =  grad_nu * (nu - 2.0_dp) * (100.0_dp - nu) / 98.0_dp

        ! normalise by n
        f = f / gt_nobs
        g = g / gt_nobs
    end subroutine garch_t_obj

    ! ---- private routines ----------------------------------------

end module garch_t_mod

! NAGARCH-X(1,1) with standardised Student-t innovations.
!
! Variance equation (same as nagarchx_mod):
!   h_t = omega + alpha*(y_{t-1} - theta*sqrt(h_{t-1}))^2 + delta*x_{t-1} + beta*h_{t-1}
!
! Innovation: y_t = sqrt(h_t) * eps_t,  eps_t ~ standardised t(nu),  nu > 2
!
! Parameterisation (6 unconstrained variables p(1:6)):
!   p(1:5)  same GARCH parameterisation as nagarchx_mod
!   p(6) -> nu = 2 + 98/(1+exp(-p(6))),  ensuring 2 < nu < 100
!
! Analytical gradient uses the same NAGARCH-X recurrences as nagarchx_mod,
! replacing the Gaussian score with the t score for h_t and adding a nu gradient:
!   factor   = (1 + score_z_std(z,dist_t,nu)*z) / (2*h_t)
!   grad_nu  = -score_shape_std(z, dist_t, nu)   [summed over t]
!   g(6)     = grad_nu * (nu-2)*(100-nu)/98       [chain rule to p(6)]

module nagarchx_t_mod
    use kind_mod,          only: dp
    use distributions_mod, only: dist_t, logpdf_std, score_z_std, score_shape_std
    use bfgs_mod,          only: bfgs_minimize
    use nagarchx_mod,      only: nagarchx_params_t, nagarchx_transform, nagarchx_inv_transform, &
                                 nagarchx_variance_path, nagarchx_persist
    implicit none
    private

    real(dp), parameter :: min_var = 1.0e-12_dp

    type, public :: nagarchx_t_params_t
        type(nagarchx_params_t) :: garch
        real(dp) :: nu = 0.0_dp
    end type nagarchx_t_params_t

    real(dp), allocatable, save :: obj_y(:), obj_x(:)
    integer,               save :: obj_ntrain = 0

    public :: fit_nagarchx_t, nagarchx_t_persist

contains

    pure real(dp) function nagarchx_t_persist(p)
        type(nagarchx_t_params_t), intent(in) :: p
        nagarchx_t_persist = nagarchx_persist(p%garch)
    end function nagarchx_t_persist

    subroutine fit_nagarchx_t(y, x, ntrain, max_iter, gtol, fopt, params, niter, converged, warm)
        real(dp),                            intent(in)  :: y(:), x(:), gtol
        integer,                             intent(in)  :: ntrain, max_iter
        real(dp),                            intent(out) :: fopt
        type(nagarchx_t_params_t),           intent(out) :: params
        integer,                             intent(out) :: niter
        logical,                             intent(out) :: converged
        type(nagarchx_t_params_t), optional, intent(in)  :: warm
        real(dp) :: pv(6), mu_y2, mu_x, nu_0
        type(nagarchx_params_t) :: garch_0

        if (size(y) /= size(x)) error stop "fit_nagarchx_t: y and x must have the same length"
        if (ntrain < 5 .or. ntrain > size(y)) error stop "fit_nagarchx_t: invalid ntrain"
        if (allocated(obj_y)) deallocate(obj_y, obj_x)
        allocate(obj_y(size(y)), obj_x(size(x)))
        obj_y      = y
        obj_x      = max(x, min_var)
        obj_ntrain = ntrain

        if (present(warm)) then
            call nagarchx_inv_transform(warm%garch, pv(1:5))
            pv(6) = log((warm%nu - 2.0_dp) / (100.0_dp - warm%nu))
        else
            nu_0  = 8.0_dp
            mu_y2 = max(sum(y(1:ntrain)**2) / real(ntrain, dp), min_var)
            mu_x  = max(sum(x(1:ntrain))    / real(ntrain, dp), min_var)
            garch_0 = nagarchx_params_t(omega = 0.05_dp * mu_y2, &
                                        alpha = 0.10_dp / 1.25_dp, &
                                        beta  = 0.80_dp, &
                                        theta = 0.5_dp, &
                                        delta = 0.05_dp * mu_y2 / mu_x)
            call nagarchx_inv_transform(garch_0, pv(1:5))
            pv(6) = log((nu_0 - 2.0_dp) / (100.0_dp - nu_0))
        end if

        call bfgs_minimize(nagarchx_t_obj, pv, 6, max_iter, gtol, fopt, niter, converged)
        call nagarchx_transform(pv(1:5), params%garch)
        params%nu = 2.0_dp + 98.0_dp / (1.0_dp + exp(-pv(6)))
    end subroutine fit_nagarchx_t

    subroutine nagarchx_t_obj(pv, np, f, g)
        integer,  intent(in)  :: np
        real(dp), intent(in)  :: pv(np)
        real(dp), intent(out) :: f, g(np)
        type(nagarchx_params_t) :: params
        real(dp) :: omega, alpha, beta, theta, delta, nu
        real(dp) :: s2, dmom, aa, D, h_unc
        real(dp) :: h, sqrth, r, kappa, z, factor
        real(dp) :: dh_dom, dh_dal, dh_dbe, dh_dth, dh_dde
        real(dp) :: grad_om, grad_al, grad_be, grad_th, grad_de, grad_nu
        integer  :: t

        call nagarchx_transform(pv(1:5), params)
        omega = params%omega;  alpha = params%alpha;  beta  = params%beta
        theta = params%theta;  delta = params%delta
        nu    = 2.0_dp + 98.0_dp / (1.0_dp + exp(-pv(6)))

        s2    = 1.0_dp + theta**2
        dmom  = 2.0_dp * theta
        aa    = alpha * s2
        D     = max(1.0_dp - aa - beta, 1.0e-8_dp)
        h_unc = omega / D

        h      = max(h_unc, min_var)
        dh_dom = 1.0_dp / D
        dh_dal = h_unc * s2 / D
        dh_dbe = h_unc / D
        dh_dth = alpha * dmom * h_unc / D
        dh_dde = 0.0_dp

        f       = 0.0_dp
        grad_om = 0.0_dp;  grad_al = 0.0_dp;  grad_be = 0.0_dp
        grad_th = 0.0_dp;  grad_de = 0.0_dp;  grad_nu = 0.0_dp

        do t = 1, obj_ntrain
            sqrth   = sqrt(max(h, min_var))
            z       = obj_y(t) / sqrth
            f       = f + 0.5_dp*log(h) - logpdf_std(z, dist_t, nu)
            factor  = (1.0_dp + score_z_std(z, dist_t, nu)*z) / (2.0_dp*h)
            grad_om = grad_om + factor * dh_dom
            grad_al = grad_al + factor * dh_dal
            grad_be = grad_be + factor * dh_dbe
            grad_th = grad_th + factor * dh_dth
            grad_de = grad_de + factor * dh_dde
            grad_nu = grad_nu - score_shape_std(z, dist_t, nu)
            ! advance recurrences
            r      = obj_y(t) - theta * sqrth
            kappa  = beta - alpha * theta * r / sqrth
            dh_dom = 1.0_dp          + kappa * dh_dom
            dh_dal = r**2            + kappa * dh_dal
            dh_dbe = h               + kappa * dh_dbe
            dh_dth = kappa * dh_dth  - 2.0_dp*alpha*r*sqrth
            dh_dde = obj_x(t)        + kappa * dh_dde
            h      = max(omega + alpha*r**2 + delta*obj_x(t) + beta*h, min_var)
        end do

        ! chain rule: constrained -> unconstrained
        g(1) =  grad_om * omega
        g(2) =  grad_al * alpha*(1.0_dp - aa) - grad_be * aa*beta
        g(3) = -grad_al * alpha*beta           + grad_be * beta*(1.0_dp - beta)
        g(4) =  grad_al * (-alpha*dmom/s2)     + grad_th
        g(5) =  grad_de * delta
        g(6) =  grad_nu * (nu - 2.0_dp) * (100.0_dp - nu) / 98.0_dp

        f = f / obj_ntrain
        g = g / obj_ntrain
    end subroutine nagarchx_t_obj

end module nagarchx_t_mod

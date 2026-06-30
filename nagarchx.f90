! NAGARCH-X(1,1) fitting module.
!
! Variance equation (in terms of constrained parameters):
!   h_t = omega + alpha*(y_{t-1} - theta*sqrt(h_{t-1}))^2 + delta*x_{t-1} + beta*h_{t-1}
!
! Parameterisation (5 unconstrained variables p(1:5)):
!   omega = exp(p1)
!   aa    = alpha*(1+theta^2) = exp(p2)/(1+exp(p2)+exp(p3))  }  aa + beta < 1
!   beta  =                     exp(p3)/(1+exp(p2)+exp(p3))  }
!   theta = p4  (unconstrained)
!   delta = exp(p5)
!
! Analytical gradient follows the NAGARCH recurrence (see nagarch.f90), adding
! one extra recurrence for delta:
!   dh_t/d_delta = x_{t-2} + kappa_{t-1}*dh_{t-1}/d_delta
! with dh_1/d_delta = 0 (h_1 does not depend on delta).

module nagarchx_mod
    use kind_mod,       only: dp
    use math_const_mod, only: log_sqrt_2pi
    use bfgs_mod,       only: bfgs_minimize
    implicit none
    private

    real(dp), parameter :: min_var = 1.0e-12_dp

    type, public :: nagarchx_params_t
        real(dp) :: omega = 0.0_dp
        real(dp) :: alpha = 0.0_dp
        real(dp) :: beta  = 0.0_dp
        real(dp) :: theta = 0.0_dp
        real(dp) :: delta = 0.0_dp
    end type nagarchx_params_t

    real(dp), allocatable, save :: obj_y(:), obj_x(:)
    integer,               save :: obj_ntrain = 0

    public :: fit_nagarchx, nagarchx_variance_path, nagarchx_persist, &
              nagarchx_transform, nagarchx_inv_transform

contains

    pure real(dp) function nagarchx_persist(p)
        type(nagarchx_params_t), intent(in) :: p
        nagarchx_persist = p%alpha*(1.0_dp + p%theta**2) + p%beta
    end function nagarchx_persist

    subroutine nagarchx_transform(pv, params)
        real(dp),                intent(in)  :: pv(5)
        type(nagarchx_params_t), intent(out) :: params
        real(dp) :: e2, e3, s, aa
        params%theta = pv(4)
        params%omega = exp(pv(1))
        params%delta = exp(pv(5))
        e2 = exp(pv(2));  e3 = exp(pv(3))
        s  = 1.0_dp + e2 + e3
        aa = e2 / s
        params%beta  = e3 / s
        params%alpha = aa / (1.0_dp + params%theta**2)
    end subroutine nagarchx_transform

    subroutine nagarchx_inv_transform(params, pv)
        type(nagarchx_params_t), intent(in)  :: params
        real(dp),                intent(out) :: pv(5)
        real(dp) :: aa, gam
        aa   = params%alpha * (1.0_dp + params%theta**2)
        gam  = max(1.0_dp - aa - params%beta, 1.0e-8_dp)
        pv(1) = log(max(params%omega, min_var))
        pv(2) = log(aa / gam)
        pv(3) = log(max(params%beta, min_var) / gam)
        pv(4) = params%theta
        pv(5) = log(max(params%delta, min_var))
    end subroutine nagarchx_inv_transform

    subroutine fit_nagarchx(y, x, ntrain, max_iter, gtol, fopt, params, niter, converged, warm)
        real(dp),                          intent(in)  :: y(:), x(:), gtol
        integer,                           intent(in)  :: ntrain, max_iter
        real(dp),                          intent(out) :: fopt
        type(nagarchx_params_t),           intent(out) :: params
        integer,                           intent(out) :: niter
        logical,                           intent(out) :: converged
        type(nagarchx_params_t), optional, intent(in)  :: warm
        real(dp) :: pv(5), mu_y2, mu_x
        type(nagarchx_params_t) :: p0

        if (size(y) /= size(x)) error stop "fit_nagarchx: y and x must have the same length"
        if (ntrain < 5 .or. ntrain > size(y)) error stop "fit_nagarchx: invalid ntrain"
        if (allocated(obj_y)) deallocate(obj_y)
        if (allocated(obj_x)) deallocate(obj_x)
        allocate(obj_y(size(y)), obj_x(size(x)))
        obj_y = y
        obj_x = max(x, min_var)
        obj_ntrain = ntrain

        if (present(warm)) then
            p0 = warm
        else
            mu_y2 = max(sum(y(1:ntrain)**2) / real(ntrain, dp), min_var)
            mu_x  = max(sum(x(1:ntrain))    / real(ntrain, dp), min_var)
            p0 = nagarchx_params_t(omega = 0.05_dp * mu_y2, &
                                   alpha = 0.10_dp / 1.25_dp, &
                                   beta  = 0.80_dp, &
                                   theta = 0.5_dp, &
                                   delta = 0.05_dp * mu_y2 / mu_x)
        end if
        call nagarchx_inv_transform(p0, pv)
        call bfgs_minimize(nagarchx_obj, pv, 5, max_iter, gtol, fopt, niter, converged)
        call nagarchx_transform(pv, params)
    end subroutine fit_nagarchx

    subroutine nagarchx_obj(pv, np, f, g)
        integer,  intent(in)  :: np
        real(dp), intent(in)  :: pv(np)
        real(dp), intent(out) :: f, g(np)
        type(nagarchx_params_t) :: params
        real(dp) :: omega, alpha, beta, theta, delta
        real(dp) :: s2, dmom, aa, D, h_unc
        real(dp) :: h, sqrth, r, kappa
        real(dp) :: dh_dom, dh_dal, dh_dbe, dh_dth, dh_dde
        real(dp) :: grad_om, grad_al, grad_be, grad_th, grad_de, factor
        integer  :: t

        call nagarchx_transform(pv, params)
        omega = params%omega;  alpha = params%alpha;  beta  = params%beta
        theta = params%theta;  delta = params%delta
        s2   = 1.0_dp + theta**2
        dmom = 2.0_dp * theta
        aa   = alpha * s2
        D    = max(1.0_dp - aa - beta, 1.0e-8_dp)
        h_unc = omega / D

        h      = max(h_unc, min_var)
        dh_dom = 1.0_dp / D
        dh_dal = h_unc * s2 / D
        dh_dbe = h_unc / D
        dh_dth = alpha * dmom * h_unc / D
        dh_dde = 0.0_dp

        f      = real(obj_ntrain, dp) * log_sqrt_2pi
        grad_om = 0.0_dp;  grad_al = 0.0_dp;  grad_be = 0.0_dp
        grad_th = 0.0_dp;  grad_de = 0.0_dp

        do t = 1, obj_ntrain
            factor  = 0.5_dp * (1.0_dp/h - obj_y(t)**2 / h**2)
            f       = f       + 0.5_dp*(log(h) + obj_y(t)**2 / h)
            grad_om = grad_om + factor * dh_dom
            grad_al = grad_al + factor * dh_dal
            grad_be = grad_be + factor * dh_dbe
            grad_th = grad_th + factor * dh_dth
            grad_de = grad_de + factor * dh_dde
            ! update recurrences: compute h_{t+1} using y_t, x_t, h_t
            sqrth  = sqrt(max(h, min_var))
            r      = obj_y(t) - theta * sqrth
            kappa  = beta - alpha * theta * r / sqrth
            dh_dom = 1.0_dp          + kappa * dh_dom
            dh_dal = r**2            + kappa * dh_dal
            dh_dbe = h               + kappa * dh_dbe
            dh_dth = kappa * dh_dth  - 2.0_dp*alpha*r*sqrth
            dh_dde = obj_x(t)        + kappa * dh_dde
            h      = max(omega + alpha*r**2 + delta*obj_x(t) + beta*h, min_var)
        end do

        ! chain rule: constrained gradients -> unconstrained gradients
        g(1) =  grad_om * omega
        g(2) =  grad_al * alpha*(1.0_dp - aa) - grad_be * aa*beta
        g(3) = -grad_al * alpha*beta           + grad_be * beta*(1.0_dp - beta)
        g(4) =  grad_al * (-alpha*dmom/s2)     + grad_th
        g(5) =  grad_de * delta

        f = f / obj_ntrain
        g = g / obj_ntrain
    end subroutine nagarchx_obj

    subroutine nagarchx_variance_path(y, x, params, h)
        real(dp),                intent(in)  :: y(:), x(:)
        type(nagarchx_params_t), intent(in)  :: params
        real(dp),                intent(out) :: h(:)
        real(dp) :: sqrth, persist
        integer  :: t

        if (size(y) /= size(x) .or. size(h) /= size(y)) &
            error stop "nagarchx_variance_path: array size mismatch"
        persist = nagarchx_persist(params)
        h(1) = max(params%omega / max(1.0_dp - persist, 1.0e-8_dp), min_var)
        do t = 2, size(y)
            sqrth = sqrt(max(h(t-1), min_var))
            h(t) = max(params%omega + params%alpha*(y(t-1) - params%theta*sqrth)**2 + &
                       params%delta*x(t-1) + params%beta*h(t-1), min_var)
        end do
    end subroutine nagarchx_variance_path

end module nagarchx_mod

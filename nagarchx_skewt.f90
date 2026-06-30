! NAGARCH-X(1,1) with Fernandez-Steel skewed-t innovations.
!
! Variance equation (same as nagarchx_mod):
!   h_t = omega + alpha*(y_{t-1} - theta*sqrt(h_{t-1}))^2 + delta*x_{t-1} + beta*h_{t-1}
!
! Innovation: y_t = sqrt(h_t) * eps_t,  eps_t ~ standardised FS-skew-t(nu, xi)
!   nu > 2 (degrees of freedom),  xi > 0 (skewness; xi=1 is symmetric)
!
! Parameterisation (7 unconstrained variables p(1:7)):
!   p(1:5)  same GARCH parameterisation as nagarchx_mod
!   p(6) -> nu = 2 + 98/(1+exp(-p(6))),  ensuring 2 < nu < 100
!   p(7) -> xi = exp(p(7)),               ensuring xi > 0
!
! Gradient is computed by central finite differences (14 NLL evaluations per BFGS step),
! consistent with how FS_SKEWT is handled elsewhere in this codebase.

module nagarchx_skewt_mod
    use kind_mod,          only: dp
    use distributions_mod, only: pdf_fs_skewt
    use bfgs_mod,          only: bfgs_minimize
    use nagarchx_mod,      only: nagarchx_params_t, nagarchx_transform, nagarchx_inv_transform, &
                                 nagarchx_persist
    implicit none
    private

    real(dp), parameter :: min_var = 1.0e-12_dp
    real(dp), parameter :: min_pdf = 1.0e-300_dp

    type, public :: nagarchx_skewt_params_t
        type(nagarchx_params_t) :: garch
        real(dp) :: nu = 0.0_dp
        real(dp) :: xi = 1.0_dp
    end type nagarchx_skewt_params_t

    real(dp), allocatable, save :: obj_y(:), obj_x(:)
    integer,               save :: obj_ntrain = 0

    public :: fit_nagarchx_skewt, nagarchx_skewt_persist

contains

    pure real(dp) function nagarchx_skewt_persist(p)
        type(nagarchx_skewt_params_t), intent(in) :: p
        nagarchx_skewt_persist = nagarchx_persist(p%garch)
    end function nagarchx_skewt_persist

    subroutine fit_nagarchx_skewt(y, x, ntrain, max_iter, gtol, fopt, params, niter, converged, warm)
        real(dp),                               intent(in)  :: y(:), x(:), gtol
        integer,                                intent(in)  :: ntrain, max_iter
        real(dp),                               intent(out) :: fopt
        type(nagarchx_skewt_params_t),          intent(out) :: params
        integer,                                intent(out) :: niter
        logical,                                intent(out) :: converged
        type(nagarchx_skewt_params_t), optional, intent(in) :: warm
        real(dp) :: pv(7), mu_y2, mu_x
        type(nagarchx_params_t) :: garch_0

        if (size(y) /= size(x)) error stop "fit_nagarchx_skewt: y and x must have the same length"
        if (ntrain < 5 .or. ntrain > size(y)) error stop "fit_nagarchx_skewt: invalid ntrain"
        if (allocated(obj_y)) deallocate(obj_y, obj_x)
        allocate(obj_y(size(y)), obj_x(size(x)))
        obj_y      = y
        obj_x      = max(x, min_var)
        obj_ntrain = ntrain

        if (present(warm)) then
            call nagarchx_inv_transform(warm%garch, pv(1:5))
            pv(6) = log((warm%nu - 2.0_dp) / (100.0_dp - warm%nu))
            pv(7) = log(max(warm%xi, min_var))
        else
            mu_y2 = max(sum(y(1:ntrain)**2) / real(ntrain, dp), min_var)
            mu_x  = max(sum(x(1:ntrain))    / real(ntrain, dp), min_var)
            garch_0 = nagarchx_params_t(omega = 0.05_dp * mu_y2, &
                                        alpha = 0.10_dp / 1.25_dp, &
                                        beta  = 0.80_dp, &
                                        theta = 0.5_dp, &
                                        delta = 0.05_dp * mu_y2 / mu_x)
            call nagarchx_inv_transform(garch_0, pv(1:5))
            pv(6) = log((8.0_dp - 2.0_dp) / (100.0_dp - 8.0_dp))  ! nu = 8
            pv(7) = 0.0_dp                                          ! xi = 1 (symmetric)
        end if

        call bfgs_minimize(nagarchx_skewt_obj, pv, 7, max_iter, gtol, fopt, niter, converged)
        call nagarchx_transform(pv(1:5), params%garch)
        params%nu = 2.0_dp + 98.0_dp / (1.0_dp + exp(-pv(6)))
        params%xi = exp(pv(7))
    end subroutine fit_nagarchx_skewt

    subroutine nagarchx_skewt_obj(pv, np, f, g)
        integer,  intent(in)  :: np
        real(dp), intent(in)  :: pv(np)
        real(dp), intent(out) :: f, g(np)
        real(dp) :: pp(np), pm(np), fp, fm, step
        integer  :: j

        f = nagarchx_skewt_nll(pv)
        do j = 1, np
            step   = 1.0e-5_dp * max(1.0_dp, abs(pv(j)))
            pp     = pv;  pp(j) = pp(j) + step
            pm     = pv;  pm(j) = pm(j) - step
            fp     = nagarchx_skewt_nll(pp)
            fm     = nagarchx_skewt_nll(pm)
            g(j)   = (fp - fm) / (2.0_dp * step)
        end do
    end subroutine nagarchx_skewt_obj

    real(dp) function nagarchx_skewt_nll(pv)
        real(dp), intent(in) :: pv(:)
        type(nagarchx_params_t) :: params
        real(dp) :: omega, alpha, beta, theta, delta, nu, xi, aa, D
        real(dp) :: h, sqrth, r, z
        integer  :: t

        call nagarchx_transform(pv(1:5), params)
        omega = params%omega;  alpha = params%alpha;  beta  = params%beta
        theta = params%theta;  delta = params%delta
        nu    = 2.0_dp + 98.0_dp / (1.0_dp + exp(-pv(6)))
        xi    = exp(pv(7))
        aa    = alpha * (1.0_dp + theta**2)
        D     = max(1.0_dp - aa - beta, 1.0e-8_dp)
        h     = max(omega / D, min_var)

        nagarchx_skewt_nll = 0.0_dp
        do t = 1, obj_ntrain
            sqrth = sqrt(max(h, min_var))
            z     = obj_y(t) / sqrth
            nagarchx_skewt_nll = nagarchx_skewt_nll + 0.5_dp*log(h) &
                                  - log(max(pdf_fs_skewt(z, nu, xi), min_pdf))
            r = obj_y(t) - theta * sqrth
            h = max(omega + alpha*r**2 + delta*obj_x(t) + beta*h, min_var)
        end do
        nagarchx_skewt_nll = nagarchx_skewt_nll / obj_ntrain
    end function nagarchx_skewt_nll

end module nagarchx_skewt_mod

! Log-linear Realized GARCH(1,1) with Normal return and measurement errors.
!
!   r_t = sqrt(h_t) z_t
!   log h_t = omega + alpha*log h_{t-1} + beta*log x_{t-1}
!   log x_t = mu + phi*log h_t + gamma1*z_t + gamma2*(z_t^2 - 1) + u_t
!
! The model is estimated by joint Gaussian quasi-likelihood. Forecast paths for
! h_t use the fitted log-h recursion and observed lagged realized variance x.
!
! Realized EGARCH (Hansen & Huang 2016) adds a leverage term to the h equation:
!   log h_t = omega + alpha*log h_{t-1} + beta*log x_{t-1} + tau*z_{t-1}
! Same measurement equation. 9 parameters vs 8. Persistence = alpha + beta*phi.

module realized_garch_mod
    use kind_mod, only: dp
    use math_const_mod, only: log_sqrt_2pi
    use bfgs_mod, only: bfgs_minimize
    implicit none
    private

    real(dp), parameter :: min_var = 1.0e-12_dp

    type, public :: realized_garch_params_t
        real(dp) :: omega = 0.0_dp
        real(dp) :: alpha = 0.10_dp
        real(dp) :: beta = 0.80_dp
        real(dp) :: mu = 0.0_dp
        real(dp) :: phi = 1.0_dp
        real(dp) :: gamma1 = 0.0_dp
        real(dp) :: gamma2 = 0.0_dp
        real(dp) :: sigma_u = 0.30_dp
    end type realized_garch_params_t

    type, public :: realized_garch_result_t
        type(realized_garch_params_t) :: params
        real(dp) :: loglik = -huge(1.0_dp)
        real(dp) :: return_loglik = -huge(1.0_dp)
        real(dp) :: meas_loglik = -huge(1.0_dp)
        real(dp) :: persist = 0.0_dp
        integer :: niter = 0
        logical :: converged = .false.
    end type realized_garch_result_t

    type, public :: realized_egarch_params_t
        real(dp) :: omega = 0.0_dp
        real(dp) :: alpha = 0.10_dp
        real(dp) :: beta = 0.80_dp
        real(dp) :: tau = 0.0_dp
        real(dp) :: mu = 0.0_dp
        real(dp) :: phi = 1.0_dp
        real(dp) :: gamma1 = 0.0_dp
        real(dp) :: gamma2 = 0.0_dp
        real(dp) :: sigma_u = 0.30_dp
    end type realized_egarch_params_t

    type, public :: realized_egarch_result_t
        type(realized_egarch_params_t) :: params
        real(dp) :: loglik = -huge(1.0_dp)
        real(dp) :: return_loglik = -huge(1.0_dp)
        real(dp) :: meas_loglik = -huge(1.0_dp)
        real(dp) :: persist = 0.0_dp
        integer :: niter = 0
        logical :: converged = .false.
    end type realized_egarch_result_t

    real(dp), allocatable, save :: obj_y(:), obj_x(:), obj_logx(:)
    integer, save :: obj_ntrain = 0
    real(dp), save :: obj_h0 = 0.0_dp

    real(dp), allocatable, save :: re_obj_y(:), re_obj_x(:), re_obj_logx(:)
    integer, save :: re_obj_ntrain = 0
    real(dp), save :: re_obj_h0 = 0.0_dp

    public :: fit_realized_garch
    public :: realized_garch_variance_path
    public :: realized_garch_return_loglik
    public :: fit_realized_egarch
    public :: realized_egarch_variance_path
    public :: realized_egarch_return_loglik

contains

    subroutine fit_realized_garch(y, x, ntrain, max_iter, gtol, result, h)
        ! Fit Realized GARCH on y(1:ntrain), x(1:ntrain), and return h for the full sample.
        real(dp), intent(in) :: y(:), x(:), gtol
        integer, intent(in) :: ntrain, max_iter
        type(realized_garch_result_t), intent(out) :: result
        real(dp), intent(out) :: h(:)
        real(dp) :: p(8), grad(8), fopt
        integer :: niter
        logical :: converged

        if (size(y) /= size(x) .or. size(h) /= size(y)) error stop "fit_realized_garch: array sizes differ"
        if (ntrain < 30 .or. ntrain > size(y)) error stop "fit_realized_garch: invalid ntrain"
        if (allocated(obj_y)) deallocate(obj_y)
        if (allocated(obj_x)) deallocate(obj_x)
        if (allocated(obj_logx)) deallocate(obj_logx)
        allocate(obj_y(size(y)), obj_x(size(x)), obj_logx(size(x)))
        obj_y = y
        obj_x = max(x, min_var)
        obj_logx = log(obj_x)
        obj_ntrain = ntrain
        obj_h0 = sum(obj_logx(1:ntrain)) / real(ntrain, dp)

        call start_params(y, obj_logx, ntrain, p)
        call bfgs_minimize(realized_garch_obj, p, 8, max_iter, gtol, fopt, niter, converged)
        call unpack_params(p, result%params)
        call realized_garch_variance_path(x, result%params, obj_h0, h)
        call realized_garch_loglik_parts(y(1:ntrain), x(1:ntrain), result%params, obj_h0, &
                                         result%return_loglik, result%meas_loglik)
        result%loglik = result%return_loglik + result%meas_loglik
        result%persist = result%params%alpha + result%params%beta*result%params%phi
        result%niter = niter
        result%converged = converged
        grad = 0.0_dp
    end subroutine fit_realized_garch

    subroutine realized_garch_variance_path(x, params, h0, h)
        ! Compute h_t for all t from the fitted log-h recursion.
        real(dp), intent(in) :: x(:), h0
        type(realized_garch_params_t), intent(in) :: params
        real(dp), intent(out) :: h(:)
        real(dp) :: logh
        integer :: t

        if (size(x) /= size(h)) error stop "realized_garch_variance_path: array sizes differ"
        logh = h0
        h(1) = max(exp(logh), min_var)
        do t = 2, size(x)
            logh = params%omega + params%alpha*logh + params%beta*log(max(x(t - 1), min_var))
            h(t) = max(exp(logh), min_var)
        end do
    end subroutine realized_garch_variance_path

    real(dp) function realized_garch_return_loglik(y, h)
        ! Gaussian return-only log likelihood for a supplied h path.
        real(dp), intent(in) :: y(:), h(:)
        integer :: t

        if (size(y) /= size(h)) error stop "realized_garch_return_loglik: array sizes differ"
        realized_garch_return_loglik = 0.0_dp
        do t = 1, size(y)
            realized_garch_return_loglik = realized_garch_return_loglik - log_sqrt_2pi - &
                                           0.5_dp*log(max(h(t), min_var)) - &
                                           0.5_dp*y(t)**2 / max(h(t), min_var)
        end do
    end function realized_garch_return_loglik

    subroutine realized_garch_obj(p, np, f, g)
        integer, intent(in) :: np
        real(dp), intent(in) :: p(np)
        real(dp), intent(out) :: f, g(np)
        real(dp), allocatable :: pp(:), pm(:)
        real(dp) :: step, fp, fm
        integer :: j

        f = realized_garch_nll(p)
        allocate(pp(np), pm(np))
        do j = 1, np
            step = 1.0e-5_dp * max(1.0_dp, abs(p(j)))
            pp = p
            pm = p
            pp(j) = pp(j) + step
            pm(j) = pm(j) - step
            fp = realized_garch_nll(pp)
            fm = realized_garch_nll(pm)
            g(j) = (fp - fm) / (2.0_dp*step)
        end do
        deallocate(pp, pm)
    end subroutine realized_garch_obj

    real(dp) function realized_garch_nll(p)
        real(dp), intent(in) :: p(:)
        type(realized_garch_params_t) :: params
        real(dp) :: ll_ret, ll_meas

        call unpack_params(p, params)
        call realized_garch_loglik_parts(obj_y(1:obj_ntrain), obj_x(1:obj_ntrain), params, obj_h0, ll_ret, ll_meas)
        realized_garch_nll = -(ll_ret + ll_meas) / real(obj_ntrain, dp)
        if (realized_garch_nll /= realized_garch_nll .or. realized_garch_nll > 1.0e29_dp) then
            realized_garch_nll = 1.0e30_dp
        end if
    end function realized_garch_nll

    subroutine realized_garch_loglik_parts(y, x, params, h0, ll_ret, ll_meas)
        real(dp), intent(in) :: y(:), x(:), h0
        type(realized_garch_params_t), intent(in) :: params
        real(dp), intent(out) :: ll_ret, ll_meas
        real(dp) :: logh, h, z, rho, u, logx
        integer :: t

        ll_ret = 0.0_dp
        ll_meas = 0.0_dp
        logh = h0
        h = max(exp(logh), min_var)
        z = y(1) / sqrt(h)
        logx = log(max(x(1), min_var))
        rho = params%gamma1*z + params%gamma2*(z**2 - 1.0_dp)
        u = logx - params%mu - params%phi*logh - rho
        ll_ret = ll_ret - log_sqrt_2pi - 0.5_dp*log(h) - 0.5_dp*z**2
        ll_meas = ll_meas - log_sqrt_2pi - log(params%sigma_u) - 0.5_dp*(u / params%sigma_u)**2
        do t = 2, size(y)
            logh = params%omega + params%alpha*logh + params%beta*log(max(x(t - 1), min_var))
            h = max(exp(logh), min_var)
            z = y(t) / sqrt(h)
            logx = log(max(x(t), min_var))
            rho = params%gamma1*z + params%gamma2*(z**2 - 1.0_dp)
            u = logx - params%mu - params%phi*logh - rho
            ll_ret = ll_ret - log_sqrt_2pi - 0.5_dp*log(h) - 0.5_dp*z**2
            ll_meas = ll_meas - log_sqrt_2pi - log(params%sigma_u) - 0.5_dp*(u / params%sigma_u)**2
        end do
    end subroutine realized_garch_loglik_parts

    subroutine start_params(y, logx, ntrain, p)
        real(dp), intent(in) :: y(:), logx(:)
        integer, intent(in) :: ntrain
        real(dp), intent(out) :: p(8)
        real(dp) :: alpha, beta, phi, logh_mean, mu, omega, sigma_u

        alpha = 0.10_dp
        beta = 0.80_dp
        phi = 1.0_dp
        logh_mean = log(max(sum(y(1:ntrain)**2) / real(ntrain, dp), min_var))
        omega = logh_mean*(1.0_dp - alpha) - beta*sum(logx(1:ntrain)) / real(ntrain, dp)
        mu = sum(logx(1:ntrain)) / real(ntrain, dp) - phi*logh_mean
        sigma_u = 0.50_dp
        call pack_params(realized_garch_params_t(omega=omega, alpha=alpha, beta=beta, mu=mu, phi=phi, &
                                                 gamma1=-0.05_dp, gamma2=0.02_dp, sigma_u=sigma_u), p)
    end subroutine start_params

    subroutine pack_params(params, p)
        type(realized_garch_params_t), intent(in) :: params
        real(dp), intent(out) :: p(8)
        real(dp) :: s, w

        p(1) = params%omega
        s = min(max(params%alpha + params%beta, 1.0e-8_dp), 0.999_dp)
        w = min(max(params%alpha / s, 1.0e-8_dp), 1.0_dp - 1.0e-8_dp)
        p(2) = log(s / (0.999_dp - s))
        p(3) = log(w / (1.0_dp - w))
        p(4) = params%mu
        p(5) = log(max(params%phi, 1.0e-8_dp))
        p(6) = params%gamma1
        p(7) = params%gamma2
        p(8) = log(max(params%sigma_u, 1.0e-8_dp))
    end subroutine pack_params

    subroutine unpack_params(p, params)
        real(dp), intent(in) :: p(:)
        type(realized_garch_params_t), intent(out) :: params
        real(dp) :: s, w

        params%omega = p(1)
        s = 0.999_dp / (1.0_dp + exp(-p(2)))
        w = 1.0_dp / (1.0_dp + exp(-p(3)))
        params%alpha = s*w
        params%beta = s*(1.0_dp - w)
        params%mu = p(4)
        params%phi = exp(p(5))
        params%gamma1 = p(6)
        params%gamma2 = p(7)
        params%sigma_u = exp(p(8))
    end subroutine unpack_params

    ! --- Realized EGARCH ---

    subroutine fit_realized_egarch(y, x, ntrain, max_iter, gtol, result, h)
        real(dp), intent(in) :: y(:), x(:), gtol
        integer, intent(in) :: ntrain, max_iter
        type(realized_egarch_result_t), intent(out) :: result
        real(dp), intent(out) :: h(:)
        real(dp) :: p(9), grad(9), fopt
        integer :: niter
        logical :: converged

        if (size(y) /= size(x) .or. size(h) /= size(y)) error stop "fit_realized_egarch: array sizes differ"
        if (ntrain < 30 .or. ntrain > size(y)) error stop "fit_realized_egarch: invalid ntrain"
        if (allocated(re_obj_y)) deallocate(re_obj_y)
        if (allocated(re_obj_x)) deallocate(re_obj_x)
        if (allocated(re_obj_logx)) deallocate(re_obj_logx)
        allocate(re_obj_y(size(y)), re_obj_x(size(x)), re_obj_logx(size(x)))
        re_obj_y = y
        re_obj_x = max(x, min_var)
        re_obj_logx = log(re_obj_x)
        re_obj_ntrain = ntrain
        re_obj_h0 = sum(re_obj_logx(1:ntrain)) / real(ntrain, dp)

        call re_start_params(y, re_obj_logx, ntrain, p)
        call bfgs_minimize(realized_egarch_obj, p, 9, max_iter, gtol, fopt, niter, converged)
        call re_unpack_params(p, result%params)
        call realized_egarch_variance_path(y, x, result%params, re_obj_h0, h)
        call realized_egarch_loglik_parts(y(1:ntrain), x(1:ntrain), result%params, re_obj_h0, &
                                          result%return_loglik, result%meas_loglik)
        result%loglik = result%return_loglik + result%meas_loglik
        result%persist = result%params%alpha + result%params%beta*result%params%phi
        result%niter = niter
        result%converged = converged
        grad = 0.0_dp
    end subroutine fit_realized_egarch

    subroutine realized_egarch_variance_path(y, x, params, h0, h)
        real(dp), intent(in) :: y(:), x(:), h0
        type(realized_egarch_params_t), intent(in) :: params
        real(dp), intent(out) :: h(:)
        real(dp) :: logh, z
        integer :: t

        if (size(y) /= size(x) .or. size(h) /= size(y)) &
            error stop "realized_egarch_variance_path: array sizes differ"
        logh = h0
        h(1) = max(exp(logh), min_var)
        z = y(1) / sqrt(h(1))
        do t = 2, size(x)
            logh = params%omega + params%alpha*logh + params%beta*log(max(x(t-1), min_var)) + params%tau*z
            h(t) = max(exp(logh), min_var)
            z = y(t) / sqrt(h(t))
        end do
    end subroutine realized_egarch_variance_path

    real(dp) function realized_egarch_return_loglik(y, h)
        real(dp), intent(in) :: y(:), h(:)
        integer :: t

        if (size(y) /= size(h)) error stop "realized_egarch_return_loglik: array sizes differ"
        realized_egarch_return_loglik = 0.0_dp
        do t = 1, size(y)
            realized_egarch_return_loglik = realized_egarch_return_loglik - log_sqrt_2pi - &
                                            0.5_dp*log(max(h(t), min_var)) - &
                                            0.5_dp*y(t)**2 / max(h(t), min_var)
        end do
    end function realized_egarch_return_loglik

    subroutine realized_egarch_obj(p, np, f, g)
        integer, intent(in) :: np
        real(dp), intent(in) :: p(np)
        real(dp), intent(out) :: f, g(np)
        real(dp), allocatable :: pp(:), pm(:)
        real(dp) :: step, fp, fm
        integer :: j

        f = realized_egarch_nll(p)
        allocate(pp(np), pm(np))
        do j = 1, np
            step = 1.0e-5_dp * max(1.0_dp, abs(p(j)))
            pp = p;  pm = p
            pp(j) = pp(j) + step
            pm(j) = pm(j) - step
            fp = realized_egarch_nll(pp)
            fm = realized_egarch_nll(pm)
            g(j) = (fp - fm) / (2.0_dp * step)
        end do
        deallocate(pp, pm)
    end subroutine realized_egarch_obj

    real(dp) function realized_egarch_nll(p)
        real(dp), intent(in) :: p(:)
        type(realized_egarch_params_t) :: params
        real(dp) :: ll_ret, ll_meas

        call re_unpack_params(p, params)
        call realized_egarch_loglik_parts(re_obj_y(1:re_obj_ntrain), re_obj_x(1:re_obj_ntrain), &
                                          params, re_obj_h0, ll_ret, ll_meas)
        realized_egarch_nll = -(ll_ret + ll_meas) / real(re_obj_ntrain, dp)
        if (realized_egarch_nll /= realized_egarch_nll .or. realized_egarch_nll > 1.0e29_dp) &
            realized_egarch_nll = 1.0e30_dp
    end function realized_egarch_nll

    subroutine realized_egarch_loglik_parts(y, x, params, h0, ll_ret, ll_meas)
        real(dp), intent(in) :: y(:), x(:), h0
        type(realized_egarch_params_t), intent(in) :: params
        real(dp), intent(out) :: ll_ret, ll_meas
        real(dp) :: logh, h, z, rho, u, logx
        integer :: t

        ll_ret = 0.0_dp
        ll_meas = 0.0_dp
        logh = h0
        h = max(exp(logh), min_var)
        z = y(1) / sqrt(h)
        logx = log(max(x(1), min_var))
        rho = params%gamma1*z + params%gamma2*(z**2 - 1.0_dp)
        u = logx - params%mu - params%phi*logh - rho
        ll_ret  = ll_ret  - log_sqrt_2pi - 0.5_dp*log(h) - 0.5_dp*z**2
        ll_meas = ll_meas - log_sqrt_2pi - log(params%sigma_u) - 0.5_dp*(u/params%sigma_u)**2
        do t = 2, size(y)
            logh = params%omega + params%alpha*logh + params%beta*log(max(x(t-1), min_var)) + params%tau*z
            h = max(exp(logh), min_var)
            z = y(t) / sqrt(h)
            logx = log(max(x(t), min_var))
            rho = params%gamma1*z + params%gamma2*(z**2 - 1.0_dp)
            u = logx - params%mu - params%phi*logh - rho
            ll_ret  = ll_ret  - log_sqrt_2pi - 0.5_dp*log(h) - 0.5_dp*z**2
            ll_meas = ll_meas - log_sqrt_2pi - log(params%sigma_u) - 0.5_dp*(u/params%sigma_u)**2
        end do
    end subroutine realized_egarch_loglik_parts

    subroutine re_start_params(y, logx, ntrain, p)
        real(dp), intent(in) :: y(:), logx(:)
        integer, intent(in) :: ntrain
        real(dp), intent(out) :: p(9)
        real(dp) :: alpha, beta, phi, logh_mean, mu, omega, sigma_u

        alpha = 0.10_dp
        beta = 0.80_dp
        phi = 1.0_dp
        logh_mean = log(max(sum(y(1:ntrain)**2) / real(ntrain, dp), min_var))
        omega = logh_mean*(1.0_dp - alpha) - beta*sum(logx(1:ntrain)) / real(ntrain, dp)
        mu = sum(logx(1:ntrain)) / real(ntrain, dp) - phi*logh_mean
        sigma_u = 0.50_dp
        call re_pack_params(realized_egarch_params_t(omega=omega, alpha=alpha, beta=beta, tau=0.0_dp, &
                                                     mu=mu, phi=phi, gamma1=-0.05_dp, gamma2=0.02_dp, &
                                                     sigma_u=sigma_u), p)
    end subroutine re_start_params

    subroutine re_pack_params(params, p)
        type(realized_egarch_params_t), intent(in) :: params
        real(dp), intent(out) :: p(9)
        real(dp) :: s, w

        p(1) = params%omega
        s = min(max(params%alpha + params%beta, 1.0e-8_dp), 0.999_dp)
        w = min(max(params%alpha / s, 1.0e-8_dp), 1.0_dp - 1.0e-8_dp)
        p(2) = log(s / (0.999_dp - s))
        p(3) = log(w / (1.0_dp - w))
        p(4) = params%tau
        p(5) = params%mu
        p(6) = log(max(params%phi, 1.0e-8_dp))
        p(7) = params%gamma1
        p(8) = params%gamma2
        p(9) = log(max(params%sigma_u, 1.0e-8_dp))
    end subroutine re_pack_params

    subroutine re_unpack_params(p, params)
        real(dp), intent(in) :: p(:)
        type(realized_egarch_params_t), intent(out) :: params
        real(dp) :: s, w

        params%omega   = p(1)
        s = 0.999_dp / (1.0_dp + exp(-p(2)))
        w = 1.0_dp / (1.0_dp + exp(-p(3)))
        params%alpha   = s*w
        params%beta    = s*(1.0_dp - w)
        params%tau     = p(4)
        params%mu      = p(5)
        params%phi     = exp(p(6))
        params%gamma1  = p(7)
        params%gamma2  = p(8)
        params%sigma_u = exp(p(9))
    end subroutine re_unpack_params

end module realized_garch_mod

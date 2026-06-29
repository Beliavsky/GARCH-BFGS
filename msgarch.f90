! 2-state Markov-Switching GARCH (Gray 1996, collapsing approximation).
!
! Each state k in {1,2} has a GARCH(1,1) variance equation:
!   h_{k,t} = omega_k + alpha_k*r_{t-1}^2 + beta_k*h_bar_{t-1}
!
! where h_bar_{t-1} = pi_{1,t-1}*h_{1,t-1} + pi_{2,t-1}*h_{2,t-1} is Gray's collapsing
! mixture (filtered-probability-weighted average of state variances).
!
! Transition matrix P = [[p11, 1-p11], [1-p22, p22]], 9 parameters total.
! dof is shared across regimes; dof > 2 enforced via dof = exp(p9) + 2.
! Persistence of state k: alpha_k + beta_k.
!
! Gradient: analytic forward-mode propagation through the filter by default.
! Toggle via set_msgarch_analytic_grad(.false.) to use central differences instead.
!
! Physical parameter ordering for analytic gradient:
!   [omega1, alpha1, beta1, omega2, alpha2, beta2, p11, p22, nu]  (9 params)

module msgarch_mod
    use kind_mod, only: dp
    use bfgs_mod, only: bfgs_minimize
    implicit none
    private

    real(dp), parameter :: min_var = 1.0e-12_dp
    real(dp), parameter :: pi_dp      = acos(-1.0_dp)
    real(dp), parameter :: log_2pi_dp = log(2.0_dp * pi_dp)

    type, public :: msgarch_params_t
        real(dp) :: omega(2) = [1.0e-6_dp, 4.0e-6_dp]
        real(dp) :: alpha(2) = [0.05_dp, 0.10_dp]
        real(dp) :: beta(2)  = [0.90_dp, 0.80_dp]
        real(dp) :: p11 = 0.97_dp
        real(dp) :: p22 = 0.97_dp
        real(dp) :: dof = 8.0_dp
    end type msgarch_params_t

    type, public :: msgarch_result_t
        type(msgarch_params_t) :: params
        real(dp) :: loglik = -huge(1.0_dp)
        real(dp) :: persist(2) = [0.0_dp, 0.0_dp]
        integer  :: niter = 0
        logical  :: converged = .false.
    end type msgarch_result_t

    real(dp), allocatable, save :: ms_obj_y(:)
    integer,  save :: ms_obj_ntrain = 0
    logical,  save :: ms_use_analytic_grad = .true.
    logical,  save :: ms_gauss_noise       = .false.

    public :: fit_msgarch, msgarch_variance_path, set_msgarch_analytic_grad, set_msgarch_gauss_noise

contains

    subroutine set_msgarch_analytic_grad(flag)
        logical, intent(in) :: flag
        ms_use_analytic_grad = flag
    end subroutine set_msgarch_analytic_grad

    subroutine set_msgarch_gauss_noise(gauss)
        logical, intent(in) :: gauss
        ms_gauss_noise = gauss
    end subroutine set_msgarch_gauss_noise

    subroutine fit_msgarch(y, ntrain, max_iter, gtol, result, h)
        real(dp), intent(in)  :: y(:), gtol
        integer,  intent(in)  :: ntrain, max_iter
        type(msgarch_result_t), intent(out) :: result
        real(dp), intent(out) :: h(:)
        integer  :: np_ms, niter
        real(dp) :: p(9), fopt
        logical  :: converged

        if (size(h) /= size(y)) error stop "fit_msgarch: h and y sizes differ"
        if (ntrain < 20 .or. ntrain > size(y)) error stop "fit_msgarch: invalid ntrain"

        if (allocated(ms_obj_y)) deallocate(ms_obj_y)
        allocate(ms_obj_y(size(y)))
        ms_obj_y      = y
        ms_obj_ntrain = ntrain

        np_ms = merge(8, 9, ms_gauss_noise)
        call ms_start_params(y(1:ntrain), p(1:np_ms))
        call bfgs_minimize(ms_obj, p(1:np_ms), np_ms, max_iter, gtol, fopt, niter, converged)
        call ms_unpack(p(1:np_ms), result%params)
        call msgarch_variance_path(y, result%params, h)
        result%loglik     = msgarch_loglik(y(1:ntrain), result%params)
        result%persist(1) = result%params%alpha(1) + result%params%beta(1)
        result%persist(2) = result%params%alpha(2) + result%params%beta(2)
        result%niter      = niter
        result%converged  = converged
    end subroutine fit_msgarch

    subroutine msgarch_variance_path(y, params, h)
        real(dp), intent(in)  :: y(:)
        type(msgarch_params_t), intent(in) :: params
        real(dp), intent(out) :: h(:)
        real(dp) :: h1, h2, h_bar, xi1, xi2, pi1, pi2, f1, f2, ft, lf1, lf2, lmax
        integer  :: t

        if (size(h) /= size(y)) error stop "msgarch_variance_path: h and y sizes differ"
        call ms_init(params, pi1, pi2, h1, h2)
        do t = 1, size(y)
            xi1  = params%p11*pi1 + (1.0_dp - params%p22)*pi2
            xi2  = 1.0_dp - xi1
            h(t) = xi1*h1 + xi2*h2
            if (ms_gauss_noise) then
                lf1 = -0.5_dp*(log_2pi_dp + log(max(h1,min_var)) + y(t)**2/max(h1,min_var))
                lf2 = -0.5_dp*(log_2pi_dp + log(max(h2,min_var)) + y(t)**2/max(h2,min_var))
            else
                lf1 = t_logdens(y(t), h1, params%dof)
                lf2 = t_logdens(y(t), h2, params%dof)
            end if
            lmax = max(lf1, lf2)
            f1   = xi1*exp(lf1 - lmax)
            f2   = xi2*exp(lf2 - lmax)
            ft   = f1 + f2
            if (ft > 0.0_dp) then
                pi1 = f1/ft
                pi2 = f2/ft
            end if
            h_bar = pi1*h1 + pi2*h2
            h1 = max(params%omega(1) + params%alpha(1)*y(t)**2 + params%beta(1)*h_bar, min_var)
            h2 = max(params%omega(2) + params%alpha(2)*y(t)**2 + params%beta(2)*h_bar, min_var)
        end do
    end subroutine msgarch_variance_path

    real(dp) function msgarch_loglik(y, params)
        real(dp), intent(in)  :: y(:)
        type(msgarch_params_t), intent(in) :: params
        real(dp) :: h1, h2, h_bar, xi1, xi2, pi1, pi2, f1, f2, ft, lf1, lf2, lmax
        integer  :: t

        call ms_init(params, pi1, pi2, h1, h2)
        msgarch_loglik = 0.0_dp
        do t = 1, size(y)
            xi1  = params%p11*pi1 + (1.0_dp - params%p22)*pi2
            xi2  = 1.0_dp - xi1
            if (ms_gauss_noise) then
                lf1 = -0.5_dp*(log_2pi_dp + log(max(h1,min_var)) + y(t)**2/max(h1,min_var))
                lf2 = -0.5_dp*(log_2pi_dp + log(max(h2,min_var)) + y(t)**2/max(h2,min_var))
            else
                lf1 = t_logdens(y(t), h1, params%dof)
                lf2 = t_logdens(y(t), h2, params%dof)
            end if
            lmax = max(lf1, lf2)
            f1   = xi1*exp(lf1 - lmax)
            f2   = xi2*exp(lf2 - lmax)
            ft   = f1 + f2
            if (ft <= 0.0_dp) then
                msgarch_loglik = -huge(1.0_dp)
                return
            end if
            msgarch_loglik = msgarch_loglik + lmax + log(ft)
            pi1   = f1/ft
            pi2   = f2/ft
            h_bar = pi1*h1 + pi2*h2
            h1 = max(params%omega(1) + params%alpha(1)*y(t)**2 + params%beta(1)*h_bar, min_var)
            h2 = max(params%omega(2) + params%alpha(2)*y(t)**2 + params%beta(2)*h_bar, min_var)
        end do
    end function msgarch_loglik

    ! Analytic forward-mode gradient of log-likelihood w.r.t. 9 physical parameters:
    !   [omega1, alpha1, beta1, omega2, alpha2, beta2, p11, p22, nu]
    ! Jacobians dh1, dh2, dpi1 (each length 9) are propagated alongside the filter state.
    ! dpi2 = -dpi1 throughout (sum-to-one constraint).
    subroutine ms_loglik_grad(y, n, params, loglik, dloglik)
        real(dp), intent(in)  :: y(:)
        integer,  intent(in)  :: n
        type(msgarch_params_t), intent(in) :: params
        real(dp), intent(out) :: loglik, dloglik(9)

        real(dp) :: h1, h2, h_bar, xi1, xi2, pi1, pi2
        real(dp) :: lf1, lf2, lmax, f1w, f2w, ft
        real(dp) :: dlf1_h, dlf2_h, dlf1_nu, dlf2_nu, u1, u2, nu_m2
        real(dp) :: dh1(9), dh2(9), dpi1(9), dhbar(9), dxi1(9), df1(9), df2(9), dft(9)
        real(dp) :: e1, e2_val
        integer  :: t

        call ms_init(params, pi1, pi2, h1, h2)
        ! Initialise initial-condition gradients (h1_0, h2_0, pi1_0 all depend on params).
        block
            real(dp) :: gam1, gam2, denom_pi
            gam1     = max(1.0_dp - params%alpha(1) - params%beta(1), 1.0e-8_dp)
            gam2     = max(1.0_dp - params%alpha(2) - params%beta(2), 1.0e-8_dp)
            denom_pi = max(2.0_dp - params%p11 - params%p22, 1.0e-10_dp)
            dh1 = 0.0_dp
            if (h1 > min_var) then
                dh1(1) = 1.0_dp / gam1              ! d(h1_0)/d(omega1)
                dh1(2) = params%omega(1) / gam1**2  ! d(h1_0)/d(alpha1)
                dh1(3) = params%omega(1) / gam1**2  ! d(h1_0)/d(beta1)
            end if
            dh2 = 0.0_dp
            if (h2 > min_var) then
                dh2(4) = 1.0_dp / gam2              ! d(h2_0)/d(omega2)
                dh2(5) = params%omega(2) / gam2**2  ! d(h2_0)/d(alpha2)
                dh2(6) = params%omega(2) / gam2**2  ! d(h2_0)/d(beta2)
            end if
            dpi1 = 0.0_dp
            dpi1(7) =  pi1 / denom_pi               ! d(pi1_0)/d(p11)
            dpi1(8) = -pi2 / denom_pi               ! d(pi1_0)/d(p22)
        end block
        loglik  = 0.0_dp;  dloglik = 0.0_dp
        nu_m2   = max(params%dof - 2.0_dp, 0.01_dp)

        do t = 1, n
            ! --- Predicted state probabilities ---
            xi1 = params%p11*pi1 + (1.0_dp - params%p22)*pi2
            xi2 = 1.0_dp - xi1

            ! dxi1 / d(physical): through dpi1 (pi2 = 1-pi1)
            ! dxi1/d(p11) += pi1, dxi1/d(p22) -= pi2
            dxi1 = (params%p11 - (1.0_dp - params%p22)) * dpi1
            dxi1(7) = dxi1(7) + pi1     ! index 7 = p11
            dxi1(8) = dxi1(8) - pi2     ! index 8 = p22

            ! --- Densities and their derivatives ---
            h1 = max(h1, min_var);  h2 = max(h2, min_var)
            if (ms_gauss_noise) then
                lf1 = -0.5_dp*(log_2pi_dp + log(h1) + y(t)**2/h1)
                lf2 = -0.5_dp*(log_2pi_dp + log(h2) + y(t)**2/h2)
                dlf1_h  = 0.5_dp/h1 * (y(t)**2/h1 - 1.0_dp)
                dlf2_h  = 0.5_dp/h2 * (y(t)**2/h2 - 1.0_dp)
                dlf1_nu = 0.0_dp;  dlf2_nu = 0.0_dp
            else
                lf1 = t_logdens(y(t), h1, params%dof)
                lf2 = t_logdens(y(t), h2, params%dof)
                u1 = y(t)**2 / (nu_m2 * h1)
                u2 = y(t)**2 / (nu_m2 * h2)
                dlf1_h  = 0.5_dp/h1 * (-1.0_dp + (params%dof + 1.0_dp)*u1/(1.0_dp + u1))
                dlf2_h  = 0.5_dp/h2 * (-1.0_dp + (params%dof + 1.0_dp)*u2/(1.0_dp + u2))
                dlf1_nu = t_logdens_dnu(y(t), h1, params%dof)
                dlf2_nu = t_logdens_dnu(y(t), h2, params%dof)
            end if

            ! --- Log-sum-exp mixture ---
            lmax = max(lf1, lf2)
            e1   = exp(lf1 - lmax);  e2_val = exp(lf2 - lmax)
            f1w  = xi1 * e1          ! unnormalised weighted density, state 1
            f2w  = xi2 * e2_val
            ft   = f1w + f2w

            if (ft <= 0.0_dp) then
                loglik = -huge(1.0_dp);  dloglik = 0.0_dp;  return
            end if
            loglik = loglik + lmax + log(ft)

            ! --- Gradient of log(ft) ---
            ! d(f1w)/d(theta_j) = e1*(dxi1(j) + xi1*(dlf1_h*dh1(j))) + e1*xi1*dlf1_nu*(j==9)
            df1 = e1 * (dxi1 + xi1*dlf1_h*dh1)
            df1(9) = df1(9) + xi1*e1*dlf1_nu          ! direct nu dependence

            ! d(f2w): dxi2 = -dxi1
            df2 = e2_val * (-dxi1 + xi2*dlf2_h*dh2)
            df2(9) = df2(9) + xi2*e2_val*dlf2_nu

            dft = df1 + df2
            dloglik = dloglik + dft/ft

            ! --- Filter update (Bayes) ---
            pi1 = f1w/ft
            pi2 = f2w/ft

            ! d(pi1)/d(theta_j) = (df1(j)*ft - f1w*dft(j)) / ft^2
            !                    = (df1(j) - pi1*dft(j)) / ft
            dpi1 = (df1 - pi1*dft) / ft

            ! --- Collapsed (filtered) mixture variance ---
            h_bar = pi1*h1 + pi2*h2
            dhbar = pi1*dh1 + pi2*dh2 + (h1 - h2)*dpi1

            ! --- Next-period state variances ---
            ! h1_new = omega1 + alpha1*y^2 + beta1*h_bar
            dh1 = params%beta(1) * dhbar
            dh1(1) = dh1(1) + 1.0_dp         ! d/d(omega1)
            dh1(2) = dh1(2) + y(t)**2         ! d/d(alpha1)
            dh1(3) = dh1(3) + h_bar           ! d/d(beta1)

            dh2 = params%beta(2) * dhbar
            dh2(4) = dh2(4) + 1.0_dp         ! d/d(omega2)
            dh2(5) = dh2(5) + y(t)**2         ! d/d(alpha2)
            dh2(6) = dh2(6) + h_bar           ! d/d(beta2)

            h1 = max(params%omega(1) + params%alpha(1)*y(t)**2 + params%beta(1)*h_bar, min_var)
            h2 = max(params%omega(2) + params%alpha(2)*y(t)**2 + params%beta(2)*h_bar, min_var)
        end do
    end subroutine ms_loglik_grad

    ! Chain rule: convert gradient w.r.t. physical params to gradient w.r.t. packed params.
    ! Physical: [omega1, alpha1, beta1, omega2, alpha2, beta2, p11, p22, nu]
    ! Packed:   [log(w1), log(a1/gam1), log(b1/gam1), log(w2), ..., logit(p11), logit(p22), log(nu-2)]
    subroutine ms_phys_to_packed_grad(params, dloglik_phys, dloglik_packed)
        type(msgarch_params_t), intent(in) :: params
        real(dp), intent(in)  :: dloglik_phys(9)
        real(dp), intent(out) :: dloglik_packed(:)  ! size 8 (Gaussian) or 9 (t)
        real(dp) :: a1, b1, a2, b2

        a1 = params%alpha(1);  b1 = params%beta(1)
        a2 = params%alpha(2);  b2 = params%beta(2)

        ! d/dp1 = d/d(omega1) * omega1
        dloglik_packed(1) = dloglik_phys(1) * params%omega(1)

        ! d/dp2 = d/d(alpha1)*alpha1*(1-alpha1) + d/d(beta1)*(-alpha1*beta1)
        ! d/dp3 = d/d(alpha1)*(-alpha1*beta1)   + d/d(beta1)*beta1*(1-beta1)
        dloglik_packed(2) = dloglik_phys(2)*a1*(1.0_dp - a1) &
                          + dloglik_phys(3)*(-a1*b1)
        dloglik_packed(3) = dloglik_phys(2)*(-a1*b1) &
                          + dloglik_phys(3)*b1*(1.0_dp - b1)

        ! d/dp4 = d/d(omega2) * omega2
        dloglik_packed(4) = dloglik_phys(4) * params%omega(2)

        dloglik_packed(5) = dloglik_phys(5)*a2*(1.0_dp - a2) &
                          + dloglik_phys(6)*(-a2*b2)
        dloglik_packed(6) = dloglik_phys(5)*(-a2*b2) &
                          + dloglik_phys(6)*b2*(1.0_dp - b2)

        ! d/dp7 = d/d(p11) * p11*(1-p11)
        dloglik_packed(7) = dloglik_phys(7) * params%p11*(1.0_dp - params%p11)
        dloglik_packed(8) = dloglik_phys(8) * params%p22*(1.0_dp - params%p22)

        if (size(dloglik_packed) >= 9) &
            dloglik_packed(9) = dloglik_phys(9) * (params%dof - 2.0_dp)
    end subroutine ms_phys_to_packed_grad

    real(dp) function t_logdens(r, h, nu)
        real(dp), intent(in) :: r, h, nu
        real(dp) :: nu_m2
        nu_m2     = nu - 2.0_dp
        t_logdens = log_gamma(0.5_dp*(nu + 1.0_dp)) - log_gamma(0.5_dp*nu) &
                   - 0.5_dp*log(nu_m2*pi_dp) - 0.5_dp*log(max(h, min_var)) &
                   - 0.5_dp*(nu + 1.0_dp)*log(1.0_dp + r**2/(nu_m2*max(h, min_var)))
    end function t_logdens

    ! d/d(nu) of log t-density at (r, h, nu).
    ! = 0.5*[digamma((nu+1)/2) - digamma(nu/2)] - 1/(2*(nu-2))
    !   - 0.5*log(1+u) + (nu+1)*u / (2*(nu-2)*(1+u))
    ! where u = r^2 / ((nu-2)*h).
    real(dp) function t_logdens_dnu(r, h, nu)
        real(dp), intent(in) :: r, h, nu
        real(dp) :: nu_m2, u
        nu_m2 = max(nu - 2.0_dp, 0.01_dp)
        u = r**2 / (nu_m2 * max(h, min_var))
        t_logdens_dnu = 0.5_dp*(digamma(0.5_dp*(nu + 1.0_dp)) - digamma(0.5_dp*nu)) &
                      - 0.5_dp/nu_m2 &
                      - 0.5_dp*log(1.0_dp + u) &
                      + 0.5_dp*(nu + 1.0_dp)/nu_m2 * u/(1.0_dp + u)
    end function t_logdens_dnu

    ! Digamma function psi(x) = d/dx log Gamma(x).
    ! Shifts x up to >= 6 via psi(x) = psi(x+1) - 1/x, then uses
    ! the asymptotic expansion: psi(z) ~ ln z - 1/(2z) - 1/(12z^2) + 1/(120z^4) - 1/(252z^6).
    real(dp) function digamma(x)
        real(dp), intent(in) :: x
        real(dp) :: z, r
        z = x;  r = 0.0_dp
        do while (z < 6.0_dp)
            r = r - 1.0_dp/z
            z = z + 1.0_dp
        end do
        r = r + log(z) - 0.5_dp/z &
              - 1.0_dp/(12.0_dp*z**2) &
              + 1.0_dp/(120.0_dp*z**4) &
              - 1.0_dp/(252.0_dp*z**6)
        digamma = r
    end function digamma

    subroutine ms_init(params, pi1, pi2, h1, h2)
        type(msgarch_params_t), intent(in) :: params
        real(dp), intent(out) :: pi1, pi2, h1, h2
        real(dp) :: denom

        denom = max(2.0_dp - params%p11 - params%p22, 1.0e-10_dp)
        pi1   = (1.0_dp - params%p22) / denom
        pi2   = 1.0_dp - pi1
        h1    = max(params%omega(1) / max(1.0_dp - params%alpha(1) - params%beta(1), 1.0e-8_dp), min_var)
        h2    = max(params%omega(2) / max(1.0_dp - params%alpha(2) - params%beta(2), 1.0e-8_dp), min_var)
    end subroutine ms_init

    subroutine ms_obj(p, np, f, g)
        integer,  intent(in)  :: np
        real(dp), intent(in)  :: p(np)
        real(dp), intent(out) :: f, g(np)
        real(dp) :: pp(np), pm(np), step, fp, fm
        real(dp) :: loglik, dloglik_phys(9)
        type(msgarch_params_t) :: params
        integer  :: j

        f = ms_nll(p)
        if (ms_use_analytic_grad) then
            call ms_unpack(p, params)
            call ms_loglik_grad(ms_obj_y, ms_obj_ntrain, params, loglik, dloglik_phys)
            call ms_phys_to_packed_grad(params, dloglik_phys, g)
            g = -g / real(ms_obj_ntrain, dp)
        else
            do j = 1, np
                step  = 1.0e-5_dp * max(1.0_dp, abs(p(j)))
                pp    = p;  pm = p
                pp(j) = pp(j) + step
                pm(j) = pm(j) - step
                fp    = ms_nll(pp)
                fm    = ms_nll(pm)
                g(j)  = (fp - fm) / (2.0_dp*step)
            end do
        end if
    end subroutine ms_obj

    real(dp) function ms_nll(p)
        real(dp), intent(in) :: p(:)
        type(msgarch_params_t) :: params

        call ms_unpack(p, params)
        ms_nll = -msgarch_loglik(ms_obj_y(1:ms_obj_ntrain), params) / real(ms_obj_ntrain, dp)
        if (ms_nll /= ms_nll .or. ms_nll > 1.0e29_dp) ms_nll = 1.0e30_dp
    end function ms_nll

    subroutine ms_pack(params, p)
        type(msgarch_params_t), intent(in)  :: params
        real(dp),               intent(out) :: p(:)  ! size 8 (Gaussian) or 9 (t)
        real(dp) :: gam

        gam  = max(1.0_dp - params%alpha(1) - params%beta(1), 1.0e-8_dp)
        p(1) = log(max(params%omega(1), min_var))
        p(2) = log(params%alpha(1) / gam)
        p(3) = log(params%beta(1)  / gam)
        gam  = max(1.0_dp - params%alpha(2) - params%beta(2), 1.0e-8_dp)
        p(4) = log(max(params%omega(2), min_var))
        p(5) = log(params%alpha(2) / gam)
        p(6) = log(params%beta(2)  / gam)
        p(7) = log(params%p11 / (1.0_dp - params%p11))
        p(8) = log(params%p22 / (1.0_dp - params%p22))
        if (size(p) >= 9) p(9) = log(max(params%dof - 2.0_dp, 0.01_dp))
    end subroutine ms_pack

    subroutine ms_unpack(p, params)
        real(dp),               intent(in)  :: p(:)  ! size 8 (Gaussian) or 9 (t)
        type(msgarch_params_t), intent(out) :: params
        real(dp) :: e2, e3, s, pv

        pv = min(max(p(1), -30.0_dp), 0.0_dp)
        params%omega(1) = exp(pv)
        e2 = exp(min(max(p(2), -20.0_dp), 20.0_dp))
        e3 = exp(min(max(p(3), -20.0_dp), 20.0_dp))
        s  = 1.0_dp + e2 + e3
        params%alpha(1) = e2/s
        params%beta(1)  = e3/s

        pv = min(max(p(4), -30.0_dp), 0.0_dp)
        params%omega(2) = exp(pv)
        e2 = exp(min(max(p(5), -20.0_dp), 20.0_dp))
        e3 = exp(min(max(p(6), -20.0_dp), 20.0_dp))
        s  = 1.0_dp + e2 + e3
        params%alpha(2) = e2/s
        params%beta(2)  = e3/s

        params%p11 = 1.0_dp / (1.0_dp + exp(-p(7)))
        params%p22 = 1.0_dp / (1.0_dp + exp(-p(8)))
        if (size(p) >= 9) then
            params%dof = exp(min(p(9), 6.0_dp)) + 2.0_dp
        else
            params%dof = 1.0e6_dp
        end if
    end subroutine ms_unpack

    subroutine ms_start_params(y, p)
        real(dp), intent(in)  :: y(:)
        real(dp), intent(out) :: p(:)  ! size 8 (Gaussian) or 9 (t)
        type(msgarch_params_t) :: params
        real(dp) :: sigma2

        sigma2 = max(sum(y**2) / real(size(y), dp), min_var)
        params%omega(1) = 0.05_dp * sigma2
        params%alpha(1) = 0.05_dp
        params%beta(1)  = 0.90_dp
        params%omega(2) = 0.20_dp * sigma2
        params%alpha(2) = 0.10_dp
        params%beta(2)  = 0.80_dp
        params%p11 = 0.97_dp
        params%p22 = 0.97_dp
        params%dof = 8.0_dp
        call ms_pack(params, p)
    end subroutine ms_start_params

end module msgarch_mod

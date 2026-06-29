! 2-state Markov-Switching NAGARCH (Gray 1996 collapsing + NAGARCH asymmetry).
!
! Each state k in {1,2} has a NAGARCH(1,1) variance equation:
!   h_{k,t} = omega_k + alpha_k*(r_{t-1} - theta_k*sqrt(h_bar_{t-1}))^2 + beta_k*h_bar_{t-1}
!
! where h_bar_{t-1} = pi_{1,t-1}*h_{1,t-1} + pi_{2,t-1}*h_{2,t-1} is Gray's collapsing
! mixture (filtered-probability-weighted average of state variances).
!
! Transition matrix P = [[p11, 1-p11], [1-p22, p22]], 11 parameters total.
! dof is shared across regimes; dof > 2 enforced via dof = exp(p11) + 2.
! Persistence of state k: alpha_k*(1 + theta_k^2) + beta_k.
!
! Gradient: analytic forward-mode propagation through the filter by default.
! Toggle via set_msnagarch_analytic_grad(.false.) to use central differences instead.
!
! Physical parameter ordering for analytic gradient:
!   [omega1, alpha1, beta1, theta1, omega2, alpha2, beta2, theta2, p11, p22, nu]  (11 params)

module msnagarch_mod
    use kind_mod, only: dp
    use bfgs_mod, only: bfgs_minimize
    implicit none
    private

    real(dp), parameter :: min_var = 1.0e-12_dp
    real(dp), parameter :: pi_dp      = acos(-1.0_dp)
    real(dp), parameter :: log_2pi_dp = log(2.0_dp * pi_dp)

    type, public :: msnagarch_params_t
        real(dp) :: omega(2) = [1.0e-6_dp, 4.0e-6_dp]
        real(dp) :: alpha(2) = [0.05_dp, 0.10_dp]
        real(dp) :: beta(2)  = [0.90_dp, 0.80_dp]
        real(dp) :: theta(2) = [0.50_dp, 0.50_dp]
        real(dp) :: p11 = 0.97_dp
        real(dp) :: p22 = 0.97_dp
        real(dp) :: dof = 8.0_dp
    end type msnagarch_params_t

    type, public :: msnagarch_result_t
        type(msnagarch_params_t) :: params
        real(dp) :: loglik = -huge(1.0_dp)
        real(dp) :: persist(2) = [0.0_dp, 0.0_dp]
        integer  :: niter = 0
        logical  :: converged = .false.
    end type msnagarch_result_t

    real(dp), allocatable, save :: msn_obj_y(:)
    integer,  save :: msn_obj_ntrain = 0
    logical,  save :: msn_use_analytic_grad = .true.
    logical,  save :: msn_gauss_noise       = .false.

    public :: fit_msnagarch, msnagarch_variance_path, set_msnagarch_analytic_grad, set_msnagarch_gauss_noise

contains

    subroutine set_msnagarch_analytic_grad(flag)
        logical, intent(in) :: flag
        msn_use_analytic_grad = flag
    end subroutine set_msnagarch_analytic_grad

    subroutine set_msnagarch_gauss_noise(gauss)
        logical, intent(in) :: gauss
        msn_gauss_noise = gauss
    end subroutine set_msnagarch_gauss_noise

    subroutine fit_msnagarch(y, ntrain, max_iter, gtol, result, h)
        real(dp), intent(in)  :: y(:), gtol
        integer,  intent(in)  :: ntrain, max_iter
        type(msnagarch_result_t), intent(out) :: result
        real(dp), intent(out) :: h(:)
        integer  :: np_msn, niter
        real(dp) :: p(11), fopt
        logical  :: converged

        if (size(h) /= size(y)) error stop "fit_msnagarch: h and y sizes differ"
        if (ntrain < 20 .or. ntrain > size(y)) error stop "fit_msnagarch: invalid ntrain"

        if (allocated(msn_obj_y)) deallocate(msn_obj_y)
        allocate(msn_obj_y(size(y)))
        msn_obj_y      = y
        msn_obj_ntrain = ntrain

        np_msn = merge(10, 11, msn_gauss_noise)
        call msn_start_params(y(1:ntrain), p(1:np_msn))
        call bfgs_minimize(msn_obj, p(1:np_msn), np_msn, max_iter, gtol, fopt, niter, converged)
        call msn_unpack(p(1:np_msn), result%params)
        call msnagarch_variance_path(y, result%params, h)
        result%loglik     = msnagarch_loglik(y(1:ntrain), result%params)
        result%persist(1) = result%params%alpha(1)*(1.0_dp + result%params%theta(1)**2) + result%params%beta(1)
        result%persist(2) = result%params%alpha(2)*(1.0_dp + result%params%theta(2)**2) + result%params%beta(2)
        result%niter      = niter
        result%converged  = converged
    end subroutine fit_msnagarch

    subroutine msnagarch_variance_path(y, params, h)
        real(dp), intent(in)  :: y(:)
        type(msnagarch_params_t), intent(in) :: params
        real(dp), intent(out) :: h(:)
        real(dp) :: h1, h2, h_bar, sqrt_hbar, xi1, xi2, pi1, pi2, f1, f2, ft, lf1, lf2, lmax, asym1, asym2
        integer  :: t

        if (size(h) /= size(y)) error stop "msnagarch_variance_path: h and y sizes differ"
        call msn_init(params, pi1, pi2, h1, h2)
        do t = 1, size(y)
            xi1  = params%p11*pi1 + (1.0_dp - params%p22)*pi2
            xi2  = 1.0_dp - xi1
            h(t) = xi1*h1 + xi2*h2
            if (msn_gauss_noise) then
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
            h_bar     = pi1*h1 + pi2*h2
            sqrt_hbar = sqrt(max(h_bar, min_var))
            asym1 = y(t) - params%theta(1)*sqrt_hbar
            asym2 = y(t) - params%theta(2)*sqrt_hbar
            h1 = max(params%omega(1) + params%alpha(1)*asym1**2 + params%beta(1)*h_bar, min_var)
            h2 = max(params%omega(2) + params%alpha(2)*asym2**2 + params%beta(2)*h_bar, min_var)
        end do
    end subroutine msnagarch_variance_path

    real(dp) function msnagarch_loglik(y, params)
        real(dp), intent(in)  :: y(:)
        type(msnagarch_params_t), intent(in) :: params
        real(dp) :: h1, h2, h_bar, sqrt_hbar, xi1, xi2, pi1, pi2, f1, f2, ft, lf1, lf2, lmax, asym1, asym2
        integer  :: t

        call msn_init(params, pi1, pi2, h1, h2)
        msnagarch_loglik = 0.0_dp
        do t = 1, size(y)
            xi1  = params%p11*pi1 + (1.0_dp - params%p22)*pi2
            xi2  = 1.0_dp - xi1
            if (msn_gauss_noise) then
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
                msnagarch_loglik = -huge(1.0_dp)
                return
            end if
            msnagarch_loglik = msnagarch_loglik + lmax + log(ft)
            pi1       = f1/ft
            pi2       = f2/ft
            h_bar     = pi1*h1 + pi2*h2
            sqrt_hbar = sqrt(max(h_bar, min_var))
            asym1 = y(t) - params%theta(1)*sqrt_hbar
            asym2 = y(t) - params%theta(2)*sqrt_hbar
            h1 = max(params%omega(1) + params%alpha(1)*asym1**2 + params%beta(1)*h_bar, min_var)
            h2 = max(params%omega(2) + params%alpha(2)*asym2**2 + params%beta(2)*h_bar, min_var)
        end do
    end function msnagarch_loglik

    ! Analytic forward-mode gradient of log-likelihood w.r.t. 11 physical parameters:
    !   [omega1, alpha1, beta1, theta1, omega2, alpha2, beta2, theta2, p11, p22, nu]
    !
    ! At each step, the h_k update depends on h_bar through:
    !   c_k = d(h_k)/d(h_bar) = -alpha_k*asym_k*theta_k/sqrt(h_bar) + beta_k
    subroutine msn_loglik_grad(y, n, params, loglik, dloglik)
        real(dp), intent(in)  :: y(:)
        integer,  intent(in)  :: n
        type(msnagarch_params_t), intent(in) :: params
        real(dp), intent(out) :: loglik, dloglik(11)

        real(dp) :: h1, h2, h_bar, sqrt_hbar, xi1, xi2, pi1, pi2
        real(dp) :: asym1, asym2, c1, c2
        real(dp) :: lf1, lf2, lmax, f1w, f2w, ft
        real(dp) :: dlf1_h, dlf2_h, dlf1_nu, dlf2_nu, u1, u2, nu_m2
        real(dp) :: dh1(11), dh2(11), dpi1(11), dhbar(11), dxi1(11), df1(11), df2(11), dft(11)
        real(dp) :: e1, e2_val
        integer  :: t

        call msn_init(params, pi1, pi2, h1, h2)
        ! Initialise initial-condition gradients.
        ! h1_0 = omega1/(1-alpha1*(1+theta1^2)-beta1), similarly for h2.
        ! pi1_0 = (1-p22)/(2-p11-p22).
        block
            real(dp) :: gam1, gam2, denom_pi
            gam1     = max(1.0_dp - params%alpha(1)*(1.0_dp + params%theta(1)**2) &
                           - params%beta(1), 1.0e-8_dp)
            gam2     = max(1.0_dp - params%alpha(2)*(1.0_dp + params%theta(2)**2) &
                           - params%beta(2), 1.0e-8_dp)
            denom_pi = max(2.0_dp - params%p11 - params%p22, 1.0e-10_dp)
            dh1 = 0.0_dp
            if (h1 > min_var) then
                dh1(1) = 1.0_dp / gam1                                               ! d(h1_0)/d(omega1)
                dh1(2) = params%omega(1)*(1.0_dp + params%theta(1)**2) / gam1**2    ! d(h1_0)/d(alpha1)
                dh1(3) = params%omega(1) / gam1**2                                   ! d(h1_0)/d(beta1)
                dh1(4) = 2.0_dp*params%omega(1)*params%alpha(1)*params%theta(1) / gam1**2 ! d(h1_0)/d(theta1)
            end if
            dh2 = 0.0_dp
            if (h2 > min_var) then
                dh2(5) = 1.0_dp / gam2                                               ! d(h2_0)/d(omega2)
                dh2(6) = params%omega(2)*(1.0_dp + params%theta(2)**2) / gam2**2    ! d(h2_0)/d(alpha2)
                dh2(7) = params%omega(2) / gam2**2                                   ! d(h2_0)/d(beta2)
                dh2(8) = 2.0_dp*params%omega(2)*params%alpha(2)*params%theta(2) / gam2**2 ! d(h2_0)/d(theta2)
            end if
            dpi1 = 0.0_dp
            dpi1(9)  =  pi1 / denom_pi   ! d(pi1_0)/d(p11)
            dpi1(10) = -pi2 / denom_pi   ! d(pi1_0)/d(p22)
        end block
        loglik  = 0.0_dp;  dloglik = 0.0_dp
        nu_m2   = max(params%dof - 2.0_dp, 0.01_dp)

        do t = 1, n
            ! --- Predicted state probabilities ---
            xi1 = params%p11*pi1 + (1.0_dp - params%p22)*pi2
            xi2 = 1.0_dp - xi1

            ! dxi1 / d(physical): indices 9=p11, 10=p22
            dxi1 = (params%p11 - (1.0_dp - params%p22)) * dpi1
            dxi1(9)  = dxi1(9)  + pi1
            dxi1(10) = dxi1(10) - pi2

            ! --- Densities and their derivatives ---
            h1 = max(h1, min_var);  h2 = max(h2, min_var)
            if (msn_gauss_noise) then
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
            f1w  = xi1 * e1
            f2w  = xi2 * e2_val
            ft   = f1w + f2w

            if (ft <= 0.0_dp) then
                loglik = -huge(1.0_dp);  dloglik = 0.0_dp;  return
            end if
            loglik = loglik + lmax + log(ft)

            ! --- Gradient of log(ft) ---
            df1 = e1 * (dxi1 + xi1*dlf1_h*dh1)
            df1(11) = df1(11) + xi1*e1*dlf1_nu       ! index 11 = nu

            df2 = e2_val * (-dxi1 + xi2*dlf2_h*dh2)
            df2(11) = df2(11) + xi2*e2_val*dlf2_nu

            dft = df1 + df2
            dloglik = dloglik + dft/ft

            ! --- Filter update (Bayes) ---
            pi1 = f1w/ft
            pi2 = f2w/ft
            dpi1 = (df1 - pi1*dft) / ft

            ! --- Collapsed (filtered) mixture variance ---
            h_bar     = pi1*h1 + pi2*h2
            sqrt_hbar = sqrt(max(h_bar, min_var))
            dhbar = pi1*dh1 + pi2*dh2 + (h1 - h2)*dpi1

            ! --- NAGARCH asymmetric shocks ---
            asym1 = y(t) - params%theta(1)*sqrt_hbar
            asym2 = y(t) - params%theta(2)*sqrt_hbar

            ! d(h_k)/d(h_bar) = -alpha_k*asym_k*theta_k/sqrt_hbar + beta_k
            c1 = -params%alpha(1)*asym1*params%theta(1)/sqrt_hbar + params%beta(1)
            c2 = -params%alpha(2)*asym2*params%theta(2)/sqrt_hbar + params%beta(2)

            ! Next-period state variances and their Jacobians
            ! h1_new = omega1 + alpha1*asym1^2 + beta1*h_bar
            dh1 = c1 * dhbar
            dh1(1) = dh1(1) + 1.0_dp                          ! d/d(omega1)
            dh1(2) = dh1(2) + asym1**2                        ! d/d(alpha1)
            dh1(3) = dh1(3) + h_bar                           ! d/d(beta1)
            dh1(4) = dh1(4) - 2.0_dp*params%alpha(1)*asym1*sqrt_hbar  ! d/d(theta1)

            dh2 = c2 * dhbar
            dh2(5) = dh2(5) + 1.0_dp                          ! d/d(omega2)
            dh2(6) = dh2(6) + asym2**2                        ! d/d(alpha2)
            dh2(7) = dh2(7) + h_bar                           ! d/d(beta2)
            dh2(8) = dh2(8) - 2.0_dp*params%alpha(2)*asym2*sqrt_hbar  ! d/d(theta2)

            h1 = max(params%omega(1) + params%alpha(1)*asym1**2 + params%beta(1)*h_bar, min_var)
            h2 = max(params%omega(2) + params%alpha(2)*asym2**2 + params%beta(2)*h_bar, min_var)
        end do
    end subroutine msn_loglik_grad

    ! Chain rule: convert gradient w.r.t. 11 physical params to 11 packed params.
    ! Physical: [omega1, alpha1, beta1, theta1, omega2, alpha2, beta2, theta2, p11, p22, nu]
    ! Packed:   [log(w1), log(as1/g1), log(b1/g1), theta1, log(w2), ..., logit(p11), logit(p22), log(nu-2)]
    ! Note: as_k = alpha_k*(1+theta_k^2) is the softmax-encoded alpha_star for state k.
    subroutine msn_phys_to_packed_grad(params, dloglik_phys, dloglik_packed)
        type(msnagarch_params_t), intent(in) :: params
        real(dp), intent(in)  :: dloglik_phys(11)
        real(dp), intent(out) :: dloglik_packed(:)  ! size 10 (Gaussian) or 11 (t)
        real(dp) :: as1, as2, b1, b2, t1sq, t2sq

        t1sq = 1.0_dp + params%theta(1)**2
        t2sq = 1.0_dp + params%theta(2)**2
        as1  = params%alpha(1) * t1sq    ! alpha_star for state 1
        as2  = params%alpha(2) * t2sq
        b1   = params%beta(1)
        b2   = params%beta(2)

        ! d/dp1 = d/d(omega1) * omega1
        dloglik_packed(1) = dloglik_phys(1) * params%omega(1)

        ! d/dp2: softmax Jacobian for (alpha_star1, beta1) then alpha1 = alpha_star1/t1sq
        !   d(alpha1)/d(p2) = d(alpha_star1)/d(p2) / t1sq = as1*(1-as1)/t1sq
        !   d(beta1)/d(p2)  = -as1*b1
        dloglik_packed(2) = dloglik_phys(2) * as1*(1.0_dp - as1)/t1sq &
                          + dloglik_phys(3) * (-as1*b1)

        ! d/dp3: d(alpha1)/d(p3) = -as1*b1/t1sq,  d(beta1)/d(p3) = b1*(1-b1)
        dloglik_packed(3) = dloglik_phys(2) * (-as1*b1/t1sq) &
                          + dloglik_phys(3) * b1*(1.0_dp - b1)

        ! d/dp4 = theta1 packed directly:
        !   d(alpha1)/d(theta1) = -2*theta1*alpha_star1/t1sq^2 = -2*theta1*alpha1/t1sq
        !   d(theta1)/d(theta1) = 1
        dloglik_packed(4) = dloglik_phys(2) * (-2.0_dp*params%theta(1)*params%alpha(1)/t1sq) &
                          + dloglik_phys(4) * 1.0_dp

        ! State 2 (indices 5-8)
        dloglik_packed(5) = dloglik_phys(5) * params%omega(2)

        dloglik_packed(6) = dloglik_phys(6) * as2*(1.0_dp - as2)/t2sq &
                          + dloglik_phys(7) * (-as2*b2)

        dloglik_packed(7) = dloglik_phys(6) * (-as2*b2/t2sq) &
                          + dloglik_phys(7) * b2*(1.0_dp - b2)

        dloglik_packed(8) = dloglik_phys(6) * (-2.0_dp*params%theta(2)*params%alpha(2)/t2sq) &
                          + dloglik_phys(8) * 1.0_dp

        ! Transition probabilities and dof
        dloglik_packed(9)  = dloglik_phys(9)  * params%p11*(1.0_dp - params%p11)
        dloglik_packed(10) = dloglik_phys(10) * params%p22*(1.0_dp - params%p22)
        if (size(dloglik_packed) >= 11) &
            dloglik_packed(11) = dloglik_phys(11) * (params%dof - 2.0_dp)
    end subroutine msn_phys_to_packed_grad

    real(dp) function t_logdens(r, h, nu)
        real(dp), intent(in) :: r, h, nu
        real(dp) :: nu_m2
        nu_m2     = nu - 2.0_dp
        t_logdens = log_gamma(0.5_dp*(nu + 1.0_dp)) - log_gamma(0.5_dp*nu) &
                   - 0.5_dp*log(nu_m2*pi_dp) - 0.5_dp*log(max(h, min_var)) &
                   - 0.5_dp*(nu + 1.0_dp)*log(1.0_dp + r**2/(nu_m2*max(h, min_var)))
    end function t_logdens

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

    subroutine msn_init(params, pi1, pi2, h1, h2)
        type(msnagarch_params_t), intent(in) :: params
        real(dp), intent(out) :: pi1, pi2, h1, h2
        real(dp) :: denom, persist1, persist2

        denom    = max(2.0_dp - params%p11 - params%p22, 1.0e-10_dp)
        pi1      = (1.0_dp - params%p22) / denom
        pi2      = 1.0_dp - pi1
        persist1 = params%alpha(1)*(1.0_dp + params%theta(1)**2) + params%beta(1)
        persist2 = params%alpha(2)*(1.0_dp + params%theta(2)**2) + params%beta(2)
        h1 = max(params%omega(1) / max(1.0_dp - persist1, 1.0e-8_dp), min_var)
        h2 = max(params%omega(2) / max(1.0_dp - persist2, 1.0e-8_dp), min_var)
    end subroutine msn_init

    subroutine msn_obj(p, np, f, g)
        integer,  intent(in)  :: np
        real(dp), intent(in)  :: p(np)
        real(dp), intent(out) :: f, g(np)
        real(dp) :: pp(np), pm(np), step, fp, fm
        real(dp) :: loglik, dloglik_phys(11)
        type(msnagarch_params_t) :: params
        integer  :: j

        f = msn_nll(p)
        if (msn_use_analytic_grad) then
            call msn_unpack(p, params)
            call msn_loglik_grad(msn_obj_y, msn_obj_ntrain, params, loglik, dloglik_phys)
            call msn_phys_to_packed_grad(params, dloglik_phys, g)
            g = -g / real(msn_obj_ntrain, dp)
        else
            do j = 1, np
                step  = 1.0e-5_dp * max(1.0_dp, abs(p(j)))
                pp    = p;  pm = p
                pp(j) = pp(j) + step
                pm(j) = pm(j) - step
                fp    = msn_nll(pp)
                fm    = msn_nll(pm)
                g(j)  = (fp - fm) / (2.0_dp*step)
            end do
        end if
    end subroutine msn_obj

    real(dp) function msn_nll(p)
        real(dp), intent(in) :: p(:)
        type(msnagarch_params_t) :: params

        call msn_unpack(p, params)
        msn_nll = -msnagarch_loglik(msn_obj_y(1:msn_obj_ntrain), params) / real(msn_obj_ntrain, dp)
        if (msn_nll /= msn_nll .or. msn_nll > 1.0e29_dp) msn_nll = 1.0e30_dp
    end function msn_nll

    subroutine msn_pack(params, p)
        type(msnagarch_params_t), intent(in)  :: params
        real(dp),                 intent(out) :: p(:)  ! size 10 (Gaussian) or 11 (t)
        real(dp) :: alpha_star, gam

        alpha_star = params%alpha(1)*(1.0_dp + params%theta(1)**2)
        gam  = max(1.0_dp - alpha_star - params%beta(1), 1.0e-8_dp)
        p(1) = log(max(params%omega(1), min_var))
        p(2) = log(alpha_star / gam)
        p(3) = log(params%beta(1) / gam)
        p(4) = params%theta(1)

        alpha_star = params%alpha(2)*(1.0_dp + params%theta(2)**2)
        gam  = max(1.0_dp - alpha_star - params%beta(2), 1.0e-8_dp)
        p(5) = log(max(params%omega(2), min_var))
        p(6) = log(alpha_star / gam)
        p(7) = log(params%beta(2) / gam)
        p(8) = params%theta(2)

        p(9)  = log(params%p11 / (1.0_dp - params%p11))
        p(10) = log(params%p22 / (1.0_dp - params%p22))
        if (size(p) >= 11) p(11) = log(max(params%dof - 2.0_dp, 0.01_dp))
    end subroutine msn_pack

    subroutine msn_unpack(p, params)
        real(dp),                 intent(in)  :: p(:)  ! size 10 (Gaussian) or 11 (t)
        type(msnagarch_params_t), intent(out) :: params
        real(dp) :: e2, e3, s, alpha_star, pv

        pv = min(max(p(1), -30.0_dp), 0.0_dp)
        params%omega(1) = exp(pv)
        e2 = exp(min(max(p(2), -20.0_dp), 20.0_dp))
        e3 = exp(min(max(p(3), -20.0_dp), 20.0_dp))
        s  = 1.0_dp + e2 + e3
        alpha_star      = e2/s
        params%beta(1)  = e3/s
        params%theta(1) = p(4)
        params%alpha(1) = alpha_star / max(1.0_dp + params%theta(1)**2, 1.0e-8_dp)

        pv = min(max(p(5), -30.0_dp), 0.0_dp)
        params%omega(2) = exp(pv)
        e2 = exp(min(max(p(6), -20.0_dp), 20.0_dp))
        e3 = exp(min(max(p(7), -20.0_dp), 20.0_dp))
        s  = 1.0_dp + e2 + e3
        alpha_star      = e2/s
        params%beta(2)  = e3/s
        params%theta(2) = p(8)
        params%alpha(2) = alpha_star / max(1.0_dp + params%theta(2)**2, 1.0e-8_dp)

        params%p11 = 1.0_dp / (1.0_dp + exp(-p(9)))
        params%p22 = 1.0_dp / (1.0_dp + exp(-p(10)))
        if (size(p) >= 11) then
            params%dof = exp(min(p(11), 6.0_dp)) + 2.0_dp
        else
            params%dof = 1.0e6_dp
        end if
    end subroutine msn_unpack

    subroutine msn_start_params(y, p)
        real(dp), intent(in)  :: y(:)
        real(dp), intent(out) :: p(:)  ! size 10 (Gaussian) or 11 (t)
        type(msnagarch_params_t) :: params
        real(dp) :: sigma2

        sigma2 = max(sum(y**2) / real(size(y), dp), min_var)
        params%omega(1) = 0.05_dp * sigma2
        params%alpha(1) = 0.05_dp
        params%beta(1)  = 0.90_dp
        params%theta(1) = 0.50_dp
        params%omega(2) = 0.20_dp * sigma2
        params%alpha(2) = 0.10_dp
        params%beta(2)  = 0.80_dp
        params%theta(2) = 0.50_dp
        params%p11 = 0.97_dp
        params%p22 = 0.97_dp
        params%dof = 8.0_dp
        call msn_pack(params, p)
    end subroutine msn_start_params

end module msnagarch_mod

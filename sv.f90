! Log-Normal Stochastic Volatility model, estimated by Quasi-Maximum Likelihood
! via the Kalman filter (Harvey, Ruiz & Shephard 1994).
!
! Model:
!   y_t   = exp(h_t/2) * eps_t,    eps_t  ~ dist(nu)  [Normal or standardised t]
!   h_t+1 = mu + phi*(h_t - mu) + sigma_eta * eta_t,  eta_t ~ N(0,1)
!   Corr(eps_t, eta_t) = rho   (leverage; rho=0 for proc_sv)
!
! QML transformation: w_t = log(y_t^2) - xi_u(dist,nu) = h_t + xi_t
!   Normal: xi_u = E[log chi^2(1)] = -1.2703628...,  Var(xi_t) = pi^2/2
!   t(nu):  xi_u(nu) = psi(1/2) + log(nu-2) - psi(nu/2)
!           Var(xi_t,nu) = pi^2/2 + trigamma(nu/2)
!   (As nu -> inf, t recovers Normal.)
!
! State-space (proc_sv):
!   measurement: w_t   = h_t + xi_t,  Var(xi_t) = var_u(dist,nu)
!   state:       h_t+1 = mu + phi*(h_t - mu) + sigma_eta*eta_t
!   initial:     h_1   ~ N(mu, sigma_eta^2/(1-phi^2))
!
! Leverage (proc_sv_lev): EKF approximation with
!   a_{t+1} = mu + phi*(a_{t|t} - mu) + rho*sigma_eta*q_t
!   P_{t+1} = phi^2*P_{t|t} + sigma_eta^2*(1-rho^2)
!   where q_t = y_t*exp(-a_{t|t}/2) approximates eps_t.
!
! Unconstrained parameters:
!   proc_sv     + dist_normal: [mu, atanh(phi), log(sigma_eta)]             np=3
!   proc_sv     + dist_t:      [mu, atanh(phi), log(sigma_eta), p_nu]       np=4
!   proc_sv_lev + dist_normal: [mu, atanh(phi), log(sigma_eta), atanh(rho)] np=4
!   proc_sv_lev + dist_t:      [mu, atanh(phi), log(sigma_eta), atanh(rho), p_nu] np=5
!   where nu = 2 + 98/(1 + exp(-p_nu)) in (2,100).
!
! Analytical gradient via Kalman filter sensitivity recursion.
! For dist_t, additional direct gradient via FD on xi_u(nu) and var_u(nu).
!
! Note: QML log-likelihoods are not directly comparable to GARCH/GAS.

module sv_mod
    use kind_mod,       only: dp
    use math_const_mod, only: pi, log_sqrt_2pi
    use random_mod,     only: random_normal, random_t_std
    use special_mod,    only: digamma, trigamma
    implicit none
    private
    public :: sv_set_data, sv_set_types, sv_np, sv_obj, sv_pred_logl, sv_skew_kurt, &
              sv_sym_inv_transform, sv_lev_inv_transform, &
              sv_t_inv_transform,   sv_lev_t_inv_transform, &
              sv_transform, sv_simulate, &
              proc_sv, proc_sv_lev, dist_normal, dist_t, &
              n_proc, n_dist, proc_names, dist_names, has_shape, model_names

    ! Moments of log(chi^2(1)) – Normal QML constants
    real(dp), parameter :: xi_u  = -1.2703628454614782_dp
    real(dp), parameter :: var_u =  4.9348022005446793_dp
    ! psi(1/2) = -gamma_EM - 2*log(2) = -1.9635100260...
    real(dp), parameter :: psi_half = -1.9635100260214235_dp

    integer, parameter :: proc_sv     = 1   ! symmetric, rho=0
    integer, parameter :: proc_sv_lev = 2   ! leverage,  rho != 0
    integer, parameter :: dist_normal = 1
    integer, parameter :: dist_t      = 2
    integer, parameter :: n_proc      = 2
    integer, parameter :: n_dist      = 2

    character(len=8),  parameter :: proc_names(n_proc) = ["SV      ", "SV-lev  "]
    character(len=8),  parameter :: dist_names(n_dist) = ["Normal  ", "t       "]
    logical,           parameter :: has_shape(n_dist)  = [.false., .true.]
    character(len=10), parameter :: model_names(n_proc, n_dist) = reshape( &
        ["SV-N      ", "SV-lev-N  ", "SV-t      ", "SV-lev-t  "], [n_proc, n_dist])

    real(dp), allocatable, save :: y_mod(:)
    integer,               save :: n_mod
    integer,               save :: iproc_mod = proc_sv
    integer,               save :: idist_mod = dist_normal

contains

    subroutine sv_set_data(y, n)
        ! Store observation data for use by sv_obj.
        real(dp), intent(in) :: y(:)
        integer,  intent(in) :: n
        if (allocated(y_mod)) deallocate(y_mod)
        allocate(y_mod(n))
        y_mod = y(1:n)
        n_mod = n
    end subroutine sv_set_data

    subroutine sv_set_types(iproc, idist)
        ! Set process and distribution type for the next sv_np / sv_obj call.
        integer, intent(in) :: iproc   ! proc_sv or proc_sv_lev
        integer, intent(in) :: idist   ! dist_normal or dist_t
        iproc_mod = iproc
        idist_mod = idist
    end subroutine sv_set_types

    function sv_np() result(np)
        ! Number of unconstrained parameters for the current proc/dist.
        integer :: np
        np = 3
        if (iproc_mod == proc_sv_lev) np = np + 1   ! atanh(rho)
        if (idist_mod == dist_t)      np = np + 1   ! p_nu
    end function sv_np

    subroutine sv_sym_inv_transform(mu, phi, sigma_eta, p)
        ! Constrained (mu, phi, sigma_eta) -> unconstrained p(3).
        real(dp), intent(in)  :: mu, phi, sigma_eta
        real(dp), intent(out) :: p(3)
        p(1) = mu
        p(2) = atanh(phi)
        p(3) = log(sigma_eta)
    end subroutine sv_sym_inv_transform

    subroutine sv_lev_inv_transform(mu, phi, sigma_eta, rho, p)
        ! Constrained (mu, phi, sigma_eta, rho) -> unconstrained p(4).
        real(dp), intent(in)  :: mu, phi, sigma_eta, rho
        real(dp), intent(out) :: p(4)
        p(1) = mu
        p(2) = atanh(phi)
        p(3) = log(sigma_eta)
        p(4) = atanh(rho)
    end subroutine sv_lev_inv_transform

    subroutine sv_t_inv_transform(mu, phi, sigma_eta, nu, p)
        ! Constrained (mu, phi, sigma_eta, nu) -> unconstrained p(4).
        real(dp), intent(in)  :: mu, phi, sigma_eta, nu
        real(dp), intent(out) :: p(4)
        p(1) = mu
        p(2) = atanh(phi)
        p(3) = log(sigma_eta)
        p(4) = log((nu-2.0_dp)/(100.0_dp-nu))   ! inv of nu = 2 + 98/(1+exp(-p))
    end subroutine sv_t_inv_transform

    subroutine sv_lev_t_inv_transform(mu, phi, sigma_eta, rho, nu, p)
        ! Constrained (mu, phi, sigma_eta, rho, nu) -> unconstrained p(5).
        real(dp), intent(in)  :: mu, phi, sigma_eta, rho, nu
        real(dp), intent(out) :: p(5)
        p(1) = mu
        p(2) = atanh(phi)
        p(3) = log(sigma_eta)
        p(4) = atanh(rho)
        p(5) = log((nu-2.0_dp)/(100.0_dp-nu))
    end subroutine sv_lev_t_inv_transform

    subroutine sv_transform(p, mu, phi, sigma_eta, rho, nu)
        ! Unconstrained p -> constrained parameters.
        ! rho=0 for proc_sv;  nu=0 for dist_normal.
        real(dp), intent(in)  :: p(:)
        real(dp), intent(out) :: mu, phi, sigma_eta, rho, nu
        integer :: i_nu
        mu        = p(1)
        phi       = tanh(p(2))
        sigma_eta = exp(p(3))
        if (iproc_mod == proc_sv_lev) then
            rho  = tanh(p(4))
            i_nu = 5
        else
            rho  = 0.0_dp
            i_nu = 4
        end if
        if (idist_mod == dist_t) then
            nu = 2.0_dp + 98.0_dp / (1.0_dp + exp(-p(i_nu)))
        else
            nu = 0.0_dp
        end if
    end subroutine sv_transform

    subroutine sv_obj(p, np, f, g)
        ! Dispatch to symmetric or leverage objective.
        integer,  intent(in)  :: np
        real(dp), intent(in)  :: p(np)
        real(dp), intent(out) :: f, g(np)
        if (iproc_mod == proc_sv_lev) then
            call sv_obj_lev(p, np, f, g)
        else
            call sv_obj_sym(p, np, f, g)
        end if
    end subroutine sv_obj

    subroutine sv_obj_sym(p, np, f, g)
        ! QML NLL and gradient for proc_sv, dist_normal or dist_t.
        integer,  intent(in)  :: np
        real(dp), intent(in)  :: p(np)
        real(dp), intent(out) :: f, g(np)

        real(dp) :: mu, phi, sigma_eta, rho_unused, nu, dnu_dpnu
        real(dp) :: xi_ut, var_ut, xi_ut_p, var_ut_p, dxi_dnu, dvar_dnu
        real(dp) :: phi2, sig2
        real(dp) :: a, P_v, F_t, v_t, K_t, L_t, e_t, r_t, w_t, a_tt
        real(dp) :: da(np), dPv(np)
        real(dp) :: logl
        integer  :: t
        real(dp), parameter :: h_fd = 1.0e-5_dp

        call sv_transform(p, mu, phi, sigma_eta, rho_unused, nu)
        phi2 = phi**2
        sig2 = sigma_eta**2

        if (idist_mod == dist_t) then
            dnu_dpnu = (nu-2.0_dp)*(100.0_dp-nu)/98.0_dp
            xi_ut    = psi_half + log(nu-2.0_dp) - digamma(0.5_dp*nu)
            var_ut   = 0.5_dp*pi**2 + trigamma(0.5_dp*nu)
            xi_ut_p  = psi_half + log(nu+h_fd-2.0_dp) - digamma(0.5_dp*(nu+h_fd))
            var_ut_p = 0.5_dp*pi**2 + trigamma(0.5_dp*(nu+h_fd))
            dxi_dnu  = (xi_ut_p  - xi_ut)  / h_fd
            dvar_dnu = (var_ut_p - var_ut) / h_fd
        else
            xi_ut    = xi_u;    var_ut   = var_u
            dxi_dnu  = 0.0_dp;  dvar_dnu = 0.0_dp;  dnu_dpnu = 0.0_dp
        end if

        a   = mu
        P_v = sig2 / (1.0_dp - phi2)

        da   = 0.0_dp;  da(1)  = 1.0_dp
        dPv  = 0.0_dp;  dPv(2) = 2.0_dp*phi*P_v;  dPv(3) = 2.0_dp*P_v

        logl = -real(n_mod, dp) * log_sqrt_2pi
        g    = 0.0_dp

        do t = 1, n_mod
            w_t = log(y_mod(t)**2) - xi_ut
            F_t = P_v + var_ut
            v_t = w_t - a
            K_t = P_v / F_t
            L_t = 1.0_dp - K_t
            e_t = v_t / F_t
            r_t = 1.0_dp/F_t - v_t**2/F_t**2

            logl = logl - 0.5_dp*(log(F_t) + v_t**2/F_t)
            g    = g + 0.5_dp*r_t*dPv - e_t*da
            if (idist_mod == dist_t) &
                g(np) = g(np) + (0.5_dp*r_t*dvar_dnu - e_t*dxi_dnu) * dnu_dpnu

            a_tt   = a + K_t*v_t
            da(2)  = (1.0_dp - phi2)*(a_tt - mu) + phi*L_t*(da(2)  + dPv(2)*v_t/F_t)
            da(3)  = phi*L_t*(da(3)  + dPv(3)*v_t/F_t)
            da(1)  = (1.0_dp - phi)  + phi*L_t*da(1)
            dPv(2) = phi2*L_t**2*dPv(2) + 2.0_dp*phi*(1.0_dp - phi2)*P_v*L_t
            dPv(3) = phi2*L_t**2*dPv(3) + 2.0_dp*sig2
            a   = mu + phi*(a_tt - mu)
            P_v = phi2*P_v*L_t + sig2
        end do

        f = -logl / n_mod
        g =     g / n_mod
    end subroutine sv_obj_sym

    subroutine sv_obj_lev(p, np, f, g)
        ! QML NLL and gradient for proc_sv_lev, dist_normal or dist_t.
        !
        ! EKF prediction:
        !   q_t      = y_t * exp(-a_{t|t}/2)
        !   a_{t+1}  = mu + phi*(a_{t|t} - mu) + rho*sigma_eta*q_t
        !   P_{t+1}  = phi^2*P_{t|t} + sigma_eta^2*(1-rho^2)
        integer,  intent(in)  :: np
        real(dp), intent(in)  :: p(np)
        real(dp), intent(out) :: f, g(np)

        real(dp) :: mu, phi, sigma_eta, rho, nu, dnu_dpnu
        real(dp) :: xi_ut, var_ut, xi_ut_p, var_ut_p, dxi_dnu, dvar_dnu
        real(dp) :: phi2, sig2, rho2, sig2_res
        real(dp) :: a, P_v, F_t, v_t, K_t, L_t, e_t, r_t, w_t, a_tt
        real(dp) :: q_t, C_t, d_tt(np)
        real(dp) :: da(np), dPv(np)
        real(dp) :: logl
        integer  :: t
        real(dp), parameter :: h_fd = 1.0e-5_dp

        call sv_transform(p, mu, phi, sigma_eta, rho, nu)
        phi2     = phi**2
        sig2     = sigma_eta**2
        rho2     = rho**2
        sig2_res = sig2 * (1.0_dp - rho2)

        if (idist_mod == dist_t) then
            dnu_dpnu = (nu-2.0_dp)*(100.0_dp-nu)/98.0_dp
            xi_ut    = psi_half + log(nu-2.0_dp) - digamma(0.5_dp*nu)
            var_ut   = 0.5_dp*pi**2 + trigamma(0.5_dp*nu)
            xi_ut_p  = psi_half + log(nu+h_fd-2.0_dp) - digamma(0.5_dp*(nu+h_fd))
            var_ut_p = 0.5_dp*pi**2 + trigamma(0.5_dp*(nu+h_fd))
            dxi_dnu  = (xi_ut_p  - xi_ut)  / h_fd
            dvar_dnu = (var_ut_p - var_ut) / h_fd
        else
            xi_ut    = xi_u;    var_ut   = var_u
            dxi_dnu  = 0.0_dp;  dvar_dnu = 0.0_dp;  dnu_dpnu = 0.0_dp
        end if

        a   = mu
        P_v = sig2 / (1.0_dp - phi2)

        da   = 0.0_dp;  da(1)  = 1.0_dp
        dPv  = 0.0_dp;  dPv(2) = 2.0_dp*phi*P_v;  dPv(3) = 2.0_dp*P_v

        logl = -real(n_mod, dp) * log_sqrt_2pi
        g    = 0.0_dp

        do t = 1, n_mod
            w_t = log(y_mod(t)**2) - xi_ut
            F_t = P_v + var_ut
            v_t = w_t - a
            K_t = P_v / F_t
            L_t = 1.0_dp - K_t
            e_t = v_t / F_t
            r_t = 1.0_dp/F_t - v_t**2/F_t**2

            logl = logl - 0.5_dp*(log(F_t) + v_t**2/F_t)
            g    = g + 0.5_dp*r_t*dPv - e_t*da
            if (idist_mod == dist_t) &
                g(np) = g(np) + (0.5_dp*r_t*dvar_dnu - e_t*dxi_dnu) * dnu_dpnu

            a_tt = a + K_t*v_t
            q_t  = y_mod(t) * exp(-0.5_dp * a_tt)
            C_t  = phi - 0.5_dp*rho*sigma_eta*q_t

            d_tt = L_t*(da + dPv*(v_t/F_t))

            da(1)  = (1.0_dp - phi)              + C_t*d_tt(1)
            da(2)  = (1.0_dp - phi2)*(a_tt - mu) + C_t*d_tt(2)
            da(3)  = rho*sigma_eta*q_t            + C_t*d_tt(3)
            da(4)  = sigma_eta*(1.0_dp-rho2)*q_t  + C_t*d_tt(4)
            dPv(1) = 0.0_dp
            dPv(2) = phi2*L_t**2*dPv(2) + 2.0_dp*phi*(1.0_dp - phi2)*P_v*L_t
            dPv(3) = phi2*L_t**2*dPv(3) + 2.0_dp*sig2_res
            dPv(4) = phi2*L_t**2*dPv(4) - 2.0_dp*sig2*rho*(1.0_dp - rho2)

            a   = mu + phi*(a_tt - mu) + rho*sigma_eta*q_t
            P_v = phi2*P_v*L_t + sig2_res
        end do

        f = -logl / n_mod
        g =     g / n_mod
    end subroutine sv_obj_lev

    subroutine sv_pred_logl(p, np, logl)
        ! Predictive log-likelihood: sum_t log p(y_t | Y_{t-1}; p).
        ! Uses KF predicted state a_t and variance P_v to form
        !   sigma2_pred = exp(a_t + 0.5*P_v)
        ! then evaluates Normal or standardised-t density for y_t.
        ! Directly comparable to GARCH log-likelihood.
        integer,  intent(in)  :: np
        real(dp), intent(in)  :: p(np)
        real(dp), intent(out) :: logl

        real(dp) :: mu, phi, sigma_eta, rho, nu
        real(dp) :: xi_ut, var_ut, lnorm_t_const
        real(dp) :: phi2, sig2, rho2, sig2_res
        real(dp) :: a, P_v, F_t, v_t, K_t, a_tt, q_t, w_t
        real(dp) :: sigma2_pred, y_t, lscore
        integer  :: t

        call sv_transform(p, mu, phi, sigma_eta, rho, nu)
        phi2 = phi**2
        sig2 = sigma_eta**2

        if (idist_mod == dist_t) then
            xi_ut  = psi_half + log(nu-2.0_dp) - digamma(0.5_dp*nu)
            var_ut = 0.5_dp*pi**2 + trigamma(0.5_dp*nu)
            lnorm_t_const = log_gamma(0.5_dp*(nu+1.0_dp)) - log_gamma(0.5_dp*nu) &
                            - 0.5_dp*log(pi*(nu-2.0_dp))
        else
            xi_ut         = xi_u
            var_ut        = var_u
            lnorm_t_const = 0.0_dp
        end if

        if (iproc_mod == proc_sv_lev) then
            rho2     = rho**2
            sig2_res = sig2 * (1.0_dp - rho2)
        else
            sig2_res = sig2
        end if

        a    = mu
        P_v  = sig2 / (1.0_dp - phi2)
        logl = 0.0_dp

        do t = 1, n_mod
            y_t = y_mod(t)

            ! predictive log-score using a_t, P_v prior to observing y_t
            sigma2_pred = exp(a + 0.5_dp*P_v)
            if (idist_mod == dist_t) then
                lscore = lnorm_t_const - 0.5_dp*log(sigma2_pred) &
                         - 0.5_dp*(nu+1.0_dp)*log(1.0_dp + y_t**2/((nu-2.0_dp)*sigma2_pred))
            else
                lscore = -log_sqrt_2pi - 0.5_dp*log(sigma2_pred) - 0.5_dp*y_t**2/sigma2_pred
            end if
            logl = logl + lscore

            ! KF update: same observation equation as sv_obj
            w_t  = log(y_t**2) - xi_ut
            F_t  = P_v + var_ut
            v_t  = w_t - a
            K_t  = P_v / F_t
            a_tt = a + K_t*v_t

            if (iproc_mod == proc_sv_lev) then
                q_t = y_t * exp(-0.5_dp*a_tt)
                a   = mu + phi*(a_tt - mu) + rho*sigma_eta*q_t
            else
                a   = mu + phi*(a_tt - mu)
            end if
            P_v = phi2*P_v*(1.0_dp - K_t) + sig2_res
        end do

    end subroutine sv_pred_logl

    subroutine sv_skew_kurt(p, np, skew, kurt)
        ! Skew and excess kurtosis of standardised predictive residuals
        ! z_t = y_t / sqrt(exp(a_t + 0.5*P_v)), using KF predicted state.
        integer,  intent(in)  :: np
        real(dp), intent(in)  :: p(np)
        real(dp), intent(out) :: skew, kurt

        real(dp), allocatable :: zz(:)
        real(dp) :: mu, phi, sigma_eta, rho, nu
        real(dp) :: xi_ut, var_ut, phi2, sig2, rho2, sig2_res
        real(dp) :: a, P_v, F_t, v_t, K_t, a_tt, q_t, w_t, sigma2_pred
        real(dp) :: zm, zv, zs, zk, dz, dz2, rn
        integer  :: t

        call sv_transform(p, mu, phi, sigma_eta, rho, nu)
        phi2 = phi**2
        sig2 = sigma_eta**2

        if (idist_mod == dist_t) then
            xi_ut  = psi_half + log(nu-2.0_dp) - digamma(0.5_dp*nu)
            var_ut = 0.5_dp*pi**2 + trigamma(0.5_dp*nu)
        else
            xi_ut  = xi_u
            var_ut = var_u
        end if

        if (iproc_mod == proc_sv_lev) then
            rho2     = rho**2
            sig2_res = sig2 * (1.0_dp - rho2)
        else
            sig2_res = sig2
        end if

        allocate(zz(n_mod))
        a   = mu
        P_v = sig2 / (1.0_dp - phi2)

        do t = 1, n_mod
            sigma2_pred = exp(a + 0.5_dp*P_v)
            zz(t) = y_mod(t) / sqrt(sigma2_pred)
            w_t  = log(y_mod(t)**2) - xi_ut
            F_t  = P_v + var_ut
            v_t  = w_t - a
            K_t  = P_v / F_t
            a_tt = a + K_t*v_t
            if (iproc_mod == proc_sv_lev) then
                q_t = y_mod(t) * exp(-0.5_dp*a_tt)
                a   = mu + phi*(a_tt - mu) + rho*sigma_eta*q_t
            else
                a   = mu + phi*(a_tt - mu)
            end if
            P_v = phi2*P_v*(1.0_dp - K_t) + sig2_res
        end do

        rn = real(n_mod, dp)
        zm = sum(zz) / rn
        zv = 0.0_dp;  zs = 0.0_dp;  zk = 0.0_dp
        do t = 1, n_mod
            dz  = zz(t) - zm
            dz2 = dz**2
            zv  = zv + dz2
            zs  = zs + dz*dz2
            zk  = zk + dz2**2
        end do
        zv   = zv / rn
        skew = (zs/rn) / zv**1.5_dp
        kurt = (zk/rn) / zv**2 - 3.0_dp

        deallocate(zz)
    end subroutine sv_skew_kurt

    subroutine sv_simulate(mu, phi, sigma_eta, rho, idist, nu, n, seed_val, y)
        ! Simulate n observations from the log-Normal SV model.
        ! h_1 drawn from stationary N(mu, sigma_eta^2/(1-phi^2)).
        ! eta_t = rho*eps_t + sqrt(1-rho^2)*zeta_t  (set rho=0 for proc_sv).
        ! eps_t is drawn from dist Normal or standardised t(nu).
        real(dp), intent(in)  :: mu, phi, sigma_eta, rho
        integer,  intent(in)  :: idist    ! dist_normal or dist_t
        real(dp), intent(in)  :: nu       ! t degrees of freedom (ignored for Normal)
        integer,  intent(in)  :: n, seed_val
        real(dp), intent(out) :: y(n)

        real(dp) :: h, eps, eta, sig_h, rho_perp
        integer  :: t, sz
        integer, allocatable :: seed_arr(:)

        call random_seed(size=sz)
        allocate(seed_arr(sz))
        seed_arr = seed_val
        call random_seed(put=seed_arr)
        deallocate(seed_arr)

        sig_h    = sigma_eta / sqrt(1.0_dp - phi**2)
        rho_perp = sqrt(1.0_dp - rho**2)
        h        = mu + sig_h * random_normal()
        do t = 1, n
            if (idist == dist_t) then
                eps = random_t_std(nu)
            else
                eps = random_normal()
            end if
            eta  = rho*eps + rho_perp*random_normal()
            y(t) = exp(0.5_dp * h) * eps
            h    = mu + phi*(h - mu) + sigma_eta*eta
        end do
    end subroutine sv_simulate

end module sv_mod

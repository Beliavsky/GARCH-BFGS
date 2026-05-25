! Flexible GARCH(1,1): any combination of 4 variance processes × 7 noise distributions.
!
! Variance processes (flx_proc):
!   proc_garch   = 1   h_t = omega + alpha*y_{t-1}^2 + beta*h_{t-1}
!   proc_nagarch = 2   h_t = omega + alpha*(y_{t-1}-theta*sqrt(h_{t-1}))^2 + beta*h_{t-1}
!   proc_gjr     = 3   h_t = omega + (alpha+gamma*I_{t-1})*y_{t-1}^2 + beta*h_{t-1}
!   proc_egarch  = 4   log(h_t) = omega + beta*log(h_{t-1}) + alpha*(|z_{t-1}|-c) + gamma*z_{t-1}
!
! Noise distributions (flx_dist):
!   dist_normal=1  dist_t=2  dist_sech=3  dist_ged=4
!   dist_laplace=5  dist_logistic=6  dist_nig=7
!
! Parameter layout: p = [p_proc(1:np_proc), p_dist(1:np_dist)]
!   np_proc: 3 for GARCH, 4 for NAGARCH/GJR/EGARCH
!   np_dist: 1 for T/GED/NIG, 0 for NORMAL/SECH/LAPLACE/LOGISTIC
!
! Key identity: NLL_t = -log f_std(z_t) + 0.5*log(h_t)
!   where z_t = y_t/sqrt(h_t), and dist_nll_score returns -log f_std(z_t).
!   For EGARCH: NLL_t = -log f_std(z_t) + 0.5*lh_t, factor = (1 + dl_dz*z)/2.
!   For others: factor = (1 + dl_dz*z)/(2*h).

module garch_flex_mod
    use kind_mod,       only: dp
    use math_const_mod, only: log_sqrt_2pi, two_pi, pi, sqrt2, sqrt3, log2, half_log2
    use garch_mod,   only: garch_transform, garch_inv_transform
    use nagarch_mod, only: nagarch_transform, nagarch_inv_transform, nagarch_set_news_impact
    use gjr_mod,     only: gjr_transform, gjr_inv_transform
    use egarch_mod,  only: egarch_transform, egarch_inv_transform
    use special_mod,    only: digamma, bessel_k01
    use distributions_mod, only: ged_lambda
    implicit none
    private

    ! Process type constants
    integer, parameter, public :: proc_garch   = 1
    integer, parameter, public :: proc_nagarch = 2
    integer, parameter, public :: proc_gjr     = 3
    integer, parameter, public :: proc_egarch  = 4
    integer, parameter, public :: proc_avgarch  = 5  ! abs-value GARCH (Taylor-Schwert): sigma_t = om+al*|y|+be*sigma
    integer, parameter, public :: proc_anagarch = 6  ! abs-value NAGARCH: sigma_t = om+al*|y-th*sigma|+be*sigma
    integer, parameter, public :: proc_agjr     = 7  ! abs-value GJR: sigma_t = om+(al+ga*I)*|y|+be*sigma
    integer, parameter, public :: proc_cnagarch = 8  ! Component NAGARCH (Engle-Lee 1999): q_t+alpha*(r²-q)+beta*(h-q)

    character(len=8), parameter, public :: proc_names(proc_cnagarch) = &
        ["GARCH   ", "NAGARCH ", "GJR     ", "EGARCH  ", &
         "AVGARCH ", "ANAGARCH", "AGJR    ", "CNAGARCH"]

    ! Distribution type constants
    integer, parameter, public :: dist_normal   = 1
    integer, parameter, public :: dist_t        = 2
    integer, parameter, public :: dist_sech     = 3
    integer, parameter, public :: dist_ged      = 4
    integer, parameter, public :: dist_laplace  = 5
    integer, parameter, public :: dist_logistic = 6
    integer, parameter, public :: dist_nig      = 7
    integer, parameter, public :: dist_skew_t   = 8

    ! Module-level saved state
    real(dp), allocatable, save :: flx_obs(:)
    integer,               save :: flx_nobs    = 0
    integer,               save :: flx_proc    = proc_garch
    integer,               save :: flx_dist    = dist_normal
    integer,               save :: flx_np_proc = 3
    integer,               save :: flx_np_dist = 0
    logical,               save :: flx_nagarch_zero_above_shift = .false.

    public :: flex_set_data, flex_set_types, flex_np, flex_obj, flex_skew_kurt, &
              flex_set_nagarch_news_impact, &
              anagarch_transform, anagarch_inv_transform, &
              cnagarch_transform, cnagarch_inv_transform

contains

    ! Store observations for use by flex_obj.
    subroutine flex_set_data(y, n)
        integer,  intent(in) :: n
        real(dp), intent(in) :: y(n)
        if (allocated(flx_obs)) deallocate(flx_obs)
        allocate(flx_obs(n))
        flx_obs  = y
        flx_nobs = n
    end subroutine flex_set_data

    ! Set process and distribution types, update np_proc and np_dist.
    subroutine flex_set_types(proc, dist)
        integer, intent(in) :: proc, dist
        flx_proc = proc
        flx_dist = dist
        ! np_proc: 3 for GARCH/AVGARCH; 4 for NAGARCH/GJR/EGARCH/ANAGARCH/AGJR; 6 for CNAGARCH
        if (proc == proc_garch .or. proc == proc_avgarch) then
            flx_np_proc = 3
        else if (proc == proc_cnagarch) then
            flx_np_proc = 6
        else
            flx_np_proc = 4
        end if
        ! np_dist: 2 for skew_t; 1 for T, GED, NIG; 0 otherwise
        if (dist == dist_skew_t) then
            flx_np_dist = 2
        else if (dist == dist_t .or. dist == dist_ged .or. dist == dist_nig) then
            flx_np_dist = 1
        else
            flx_np_dist = 0
        end if
    end subroutine flex_set_types

    subroutine flex_set_nagarch_news_impact(zero_above_shift)
        logical, intent(in) :: zero_above_shift
        flx_nagarch_zero_above_shift = zero_above_shift
        call nagarch_set_news_impact(zero_above_shift)
    end subroutine flex_set_nagarch_news_impact

    ! Return total number of parameters.
    integer function flex_np()
        flex_np = flx_np_proc + flx_np_dist
    end function flex_np

    ! NLL/n and gradient w.r.t. unconstrained parameters.
    subroutine flex_obj(p, np, f, g)
        integer,  intent(in)  :: np
        real(dp), intent(in)  :: p(np)
        real(dp), intent(out) :: f
        real(dp), intent(out) :: g(np)

        real(dp) :: grad_proc(6), grad_dist1, grad_dist2
        real(dp) :: p_dist1, p_dist2
        real(dp) :: nll_t, dl_dz, dp1, dp2, factor, z
        integer  :: npp, nd, t

        ! Initialise accumulators
        f           = 0.0_dp
        grad_proc   = 0.0_dp
        grad_dist1  = 0.0_dp
        grad_dist2  = 0.0_dp

        ! Extract sizes and distribution parameters
        npp     = flx_np_proc
        nd      = flx_np_dist
        p_dist1 = 0.0_dp;  if (nd > 0) p_dist1 = p(npp+1)
        p_dist2 = 0.0_dp;  if (nd > 1) p_dist2 = p(npp+2)

        select case (flx_proc)

        ! ----------------------------------------------------------------
        case (proc_garch)
            block
                real(dp) :: omega, alpha, beta
                real(dp) :: h, dh_dom, dh_dal, dh_dbe
                real(dp) :: h_unc, D

                call garch_transform(p(1:3), omega, alpha, beta)
                D     = 1.0_dp - alpha - beta
                h_unc = omega / D

                h      = h_unc
                dh_dom = 1.0_dp / D
                dh_dal = h_unc / D
                dh_dbe = h_unc / D

                do t = 1, flx_nobs
                    z = flx_obs(t) / sqrt(h)
                    call dist_nll_score(flx_dist, z, p_dist1, p_dist2, nll_t, dl_dz, dp1, dp2)
                    factor = (1.0_dp + dl_dz * z) / (2.0_dp * h)
                    f           = f           + nll_t + 0.5_dp * log(h)
                    grad_proc(1) = grad_proc(1) + factor * dh_dom
                    grad_proc(2) = grad_proc(2) + factor * dh_dal
                    grad_proc(3) = grad_proc(3) + factor * dh_dbe
                    grad_dist1   = grad_dist1   + dp1
                    grad_dist2   = grad_dist2   + dp2

                    ! Advance recurrences
                    dh_dom = 1.0_dp               + beta * dh_dom
                    dh_dal = flx_obs(t)**2         + beta * dh_dal
                    dh_dbe = h                     + beta * dh_dbe
                    h      = omega + alpha * flx_obs(t)**2 + beta * h
                end do

                ! Chain rule: softmax Jacobian for (omega, alpha, beta)
                g(1) =  grad_proc(1) * omega
                g(2) =  grad_proc(2) * alpha*(1.0_dp - alpha) &
                      - grad_proc(3) * alpha*beta
                g(3) = -grad_proc(2) * alpha*beta &
                      + grad_proc(3) * beta*(1.0_dp - beta)
            end block

        ! ----------------------------------------------------------------
        case (proc_nagarch)
            block
                real(dp) :: omega, alpha, beta, theta
                real(dp) :: h, dh_dom, dh_dal, dh_dbe, dh_dth
                real(dp) :: h_unc, D, s2, aa, dmom
                real(dp) :: sqrth, r, kappa
                logical  :: active

                call nagarch_transform(p(1:4), omega, alpha, beta, theta)
                call nagarch_shift_moments(theta, s2, dmom)
                D     = 1.0_dp - alpha*s2 - beta
                h_unc = omega / D

                h      = h_unc
                dh_dom = 1.0_dp / D
                dh_dal = h_unc * s2 / D
                dh_dbe = h_unc / D
                dh_dth = alpha * dmom * h_unc / D

                do t = 1, flx_nobs
                    z = flx_obs(t) / sqrt(h)
                    call dist_nll_score(flx_dist, z, p_dist1, p_dist2, nll_t, dl_dz, dp1, dp2)
                    factor = (1.0_dp + dl_dz * z) / (2.0_dp * h)
                    f           = f           + nll_t + 0.5_dp * log(h)
                    grad_proc(1) = grad_proc(1) + factor * dh_dom
                    grad_proc(2) = grad_proc(2) + factor * dh_dal
                    grad_proc(3) = grad_proc(3) + factor * dh_dbe
                    grad_proc(4) = grad_proc(4) + factor * dh_dth
                    grad_dist1   = grad_dist1   + dp1
                    grad_dist2   = grad_dist2   + dp2

                    ! Advance recurrences
                    sqrth  = sqrt(h)
                    r      = flx_obs(t) - theta * sqrth
                    active = .true.
                    if (flx_nagarch_zero_above_shift .and. r > 0.0_dp) then
                        r = 0.0_dp
                        active = .false.
                    end if
                    kappa  = beta
                    if (active) kappa = beta - alpha * theta * r / sqrth
                    dh_dom = 1.0_dp                  + kappa * dh_dom
                    dh_dal = r**2                    + kappa * dh_dal
                    dh_dbe = h                       + kappa * dh_dbe
                    dh_dth = kappa * dh_dth
                    if (active) dh_dth = dh_dth - 2.0_dp*alpha*r*sqrth
                    h      = omega + alpha*r**2 + beta*h
                end do

                ! Chain rule for NAGARCH
                aa   = alpha * s2
                g(1) =  grad_proc(1) * omega
                g(2) =  grad_proc(2) * alpha*(1.0_dp - aa) &
                      - grad_proc(3) * aa*beta
                g(3) = -grad_proc(2) * alpha*beta &
                      + grad_proc(3) * beta*(1.0_dp - beta)
                g(4) =  grad_proc(2) * (-alpha*dmom/s2) + grad_proc(4)
            end block

        ! ----------------------------------------------------------------
        case (proc_gjr)
            block
                real(dp) :: omega, alpha, gamma, beta
                real(dp) :: h, dh_dom, dh_dal, dh_dga, dh_dbe
                real(dp) :: h_unc, D, g2, ind

                call gjr_transform(p(1:4), omega, alpha, gamma, beta)
                g2    = 0.5_dp * gamma
                D     = 1.0_dp - alpha - g2 - beta
                h_unc = omega / D

                h      = h_unc
                dh_dom = 1.0_dp / D
                dh_dal = h_unc / D
                dh_dga = h_unc / (2.0_dp * D)
                dh_dbe = h_unc / D

                do t = 1, flx_nobs
                    z = flx_obs(t) / sqrt(h)
                    call dist_nll_score(flx_dist, z, p_dist1, p_dist2, nll_t, dl_dz, dp1, dp2)
                    factor = (1.0_dp + dl_dz * z) / (2.0_dp * h)
                    f           = f           + nll_t + 0.5_dp * log(h)
                    grad_proc(1) = grad_proc(1) + factor * dh_dom
                    grad_proc(2) = grad_proc(2) + factor * dh_dal
                    grad_proc(3) = grad_proc(3) + factor * dh_dga
                    grad_proc(4) = grad_proc(4) + factor * dh_dbe
                    grad_dist1   = grad_dist1   + dp1
                    grad_dist2   = grad_dist2   + dp2

                    ! Advance recurrences
                    if (flx_obs(t) < 0.0_dp) then
                        ind = 1.0_dp
                    else
                        ind = 0.0_dp
                    end if
                    dh_dom = 1.0_dp                       + beta * dh_dom
                    dh_dal = flx_obs(t)**2                 + beta * dh_dal
                    dh_dga = ind * flx_obs(t)**2           + beta * dh_dga
                    dh_dbe = h                             + beta * dh_dbe
                    h      = omega + (alpha + gamma*ind)*flx_obs(t)**2 + beta*h
                end do

                ! Chain rule: 3-way softmax for (alpha, gamma/2, beta)
                g(1) =  grad_proc(1) * omega
                g(2) =  alpha * ( grad_proc(2)*(1.0_dp - alpha) &
                                - 2.0_dp*grad_proc(3)*g2 &
                                - grad_proc(4)*beta )
                g(3) =  g2    * ( -grad_proc(2)*alpha &
                                + 2.0_dp*grad_proc(3)*(1.0_dp - g2) &
                                - grad_proc(4)*beta )
                g(4) =  beta  * ( -grad_proc(2)*alpha &
                                - 2.0_dp*grad_proc(3)*g2 &
                                + grad_proc(4)*(1.0_dp - beta) )
            end block

        ! ----------------------------------------------------------------
        case (proc_egarch)
            block
                real(dp) :: omega, alpha, gamma, beta
                real(dp) :: lh, h, abs_z, kappa, c_eg
                real(dp) :: dlh_dom, dlh_dal, dlh_dga, dlh_dbe

                call egarch_transform(p(1:4), omega, alpha, gamma, beta)

                ! c_eg = sqrt(2/pi), precomputed once before the loop
                c_eg = sqrt(2.0_dp / pi)

                lh = omega / (1.0_dp - beta)

                ! Initial log-variance derivatives from lh_1 = omega/(1-beta)
                dlh_dom = 1.0_dp / (1.0_dp - beta)
                dlh_dal = 0.0_dp
                dlh_dga = 0.0_dp
                dlh_dbe = lh / (1.0_dp - beta)

                do t = 1, flx_nobs
                    h = exp(lh)
                    z = flx_obs(t) / sqrt(h)
                    call dist_nll_score(flx_dist, z, p_dist1, p_dist2, nll_t, dl_dz, dp1, dp2)
                    ! For EGARCH: factor has no division by h
                    factor = (1.0_dp + dl_dz * z) * 0.5_dp
                    f           = f           + nll_t + 0.5_dp * lh
                    grad_proc(1) = grad_proc(1) + factor * dlh_dom
                    grad_proc(2) = grad_proc(2) + factor * dlh_dal
                    grad_proc(3) = grad_proc(3) + factor * dlh_dga
                    grad_proc(4) = grad_proc(4) + factor * dlh_dbe
                    grad_dist1   = grad_dist1   + dp1
                    grad_dist2   = grad_dist2   + dp2

                    ! Advance recurrences
                    abs_z  = abs(z)
                    kappa  = beta - 0.5_dp * (alpha * abs_z + gamma * z)
                    dlh_dom = 1.0_dp               + kappa * dlh_dom
                    dlh_dal = (abs_z - c_eg)       + kappa * dlh_dal
                    dlh_dga = z                    + kappa * dlh_dga
                    dlh_dbe = lh                   + kappa * dlh_dbe
                    lh      = omega + beta*lh + alpha*(abs_z - c_eg) + gamma*z
                end do

                ! Chain rule: omega, alpha, gamma free; beta = tanh(p4)
                g(1) = grad_proc(1)
                g(2) = grad_proc(2)
                g(3) = grad_proc(3)
                g(4) = grad_proc(4) * (1.0_dp - beta**2)
            end block

        ! ----------------------------------------------------------------
        ! Abs-value models: update sigma_t directly (not variance h_t).
        ! NLL_t = nll_z(z_t) + log(sigma_t),  z_t = y_t/sigma_t
        ! factor = (1 + dl_dz*z) / sigma_t  (vs (1+dl_dz*z)/(2*h) for variance models)
        ! ----------------------------------------------------------------
        case (proc_avgarch)
            ! sigma_t = omega + alpha*|y_{t-1}| + beta*sigma_{t-1}
            ! Same 3-param softmax as GARCH: alpha+beta < 1 ensures stationarity.
            block
                real(dp) :: omega, alpha, beta
                real(dp) :: sig, ds_dom, ds_dal, ds_dbe
                real(dp) :: sig_unc, D

                call garch_transform(p(1:3), omega, alpha, beta)
                D       = 1.0_dp - alpha - beta
                sig_unc = omega / D

                sig    = sig_unc
                ds_dom = 1.0_dp / D
                ds_dal = sig_unc / D
                ds_dbe = sig_unc / D

                do t = 1, flx_nobs
                    z = flx_obs(t) / sig
                    call dist_nll_score(flx_dist, z, p_dist1, p_dist2, nll_t, dl_dz, dp1, dp2)
                    factor = (1.0_dp + dl_dz * z) / sig
                    f            = f            + nll_t + log(sig)
                    grad_proc(1) = grad_proc(1) + factor * ds_dom
                    grad_proc(2) = grad_proc(2) + factor * ds_dal
                    grad_proc(3) = grad_proc(3) + factor * ds_dbe
                    grad_dist1   = grad_dist1   + dp1
                    grad_dist2   = grad_dist2   + dp2

                    ! Advance recurrences: sigma_{t+1} = omega + alpha*|y_t| + beta*sigma_t
                    ds_dom = 1.0_dp              + beta * ds_dom
                    ds_dal = abs(flx_obs(t))     + beta * ds_dal
                    ds_dbe = sig                 + beta * ds_dbe
                    sig    = omega + alpha * abs(flx_obs(t)) + beta * sig
                end do

                ! Chain rule: same softmax Jacobian as GARCH
                g(1) =  grad_proc(1) * omega
                g(2) =  grad_proc(2) * alpha*(1.0_dp - alpha) &
                      - grad_proc(3) * alpha*beta
                g(3) = -grad_proc(2) * alpha*beta &
                      + grad_proc(3) * beta*(1.0_dp - beta)
            end block

        ! ----------------------------------------------------------------
        case (proc_anagarch)
            ! sigma_t = omega + alpha*|y_{t-1} - theta*sigma_{t-1}| + beta*sigma_{t-1}
            ! 3-param softmax on (alpha, beta); theta free.
            ! kappa_t = beta - alpha*theta*sign(r_t)  (multiplier on dsigma/dp recurrence)
            block
                real(dp) :: omega, alpha, beta, theta
                real(dp) :: sig, ds_dom, ds_dal, ds_dbe, ds_dth
                real(dp) :: sig_unc, D, r, s_r, kappa

                call anagarch_transform(p(1:4), omega, alpha, beta, theta)
                D       = 1.0_dp - alpha - beta
                sig_unc = omega / D

                sig    = sig_unc
                ds_dom = 1.0_dp / D
                ds_dal = sig_unc / D
                ds_dbe = sig_unc / D
                ds_dth = 0.0_dp

                do t = 1, flx_nobs
                    z = flx_obs(t) / sig
                    call dist_nll_score(flx_dist, z, p_dist1, p_dist2, nll_t, dl_dz, dp1, dp2)
                    factor = (1.0_dp + dl_dz * z) / sig
                    f            = f            + nll_t + log(sig)
                    grad_proc(1) = grad_proc(1) + factor * ds_dom
                    grad_proc(2) = grad_proc(2) + factor * ds_dal
                    grad_proc(3) = grad_proc(3) + factor * ds_dbe
                    grad_proc(4) = grad_proc(4) + factor * ds_dth
                    grad_dist1   = grad_dist1   + dp1
                    grad_dist2   = grad_dist2   + dp2

                    ! Advance recurrences: r_t = y_t - theta*sigma_t; s_r = sign(r_t)
                    r     = flx_obs(t) - theta * sig
                    s_r   = sign(1.0_dp, r)
                    kappa = beta - alpha * theta * s_r
                    ds_dom = 1.0_dp          + kappa * ds_dom
                    ds_dal = abs(r)          + kappa * ds_dal
                    ds_dbe = sig             + kappa * ds_dbe
                    ds_dth = -alpha*s_r*sig  + kappa * ds_dth
                    sig    = omega + alpha * abs(r) + beta * sig
                end do

                ! Chain rule: simple softmax for (alpha, beta); theta free (p4)
                g(1) =  grad_proc(1) * omega
                g(2) =  grad_proc(2) * alpha*(1.0_dp - alpha) &
                      - grad_proc(3) * alpha*beta
                g(3) = -grad_proc(2) * alpha*beta &
                      + grad_proc(3) * beta*(1.0_dp - beta)
                g(4) =  grad_proc(4)
            end block

        ! ----------------------------------------------------------------
        case (proc_agjr)
            ! sigma_t = omega + (alpha + gamma*I_{t-1})*|y_{t-1}| + beta*sigma_{t-1}
            ! Same 4-param GJR softmax: alpha + gamma/2 + beta < 1.
            block
                real(dp) :: omega, alpha, gamma, beta
                real(dp) :: sig, ds_dom, ds_dal, ds_dga, ds_dbe
                real(dp) :: sig_unc, D, g2, ind

                call gjr_transform(p(1:4), omega, alpha, gamma, beta)
                g2      = 0.5_dp * gamma
                D       = 1.0_dp - alpha - g2 - beta
                sig_unc = omega / D

                sig    = sig_unc
                ds_dom = 1.0_dp / D
                ds_dal = sig_unc / D
                ds_dga = sig_unc / (2.0_dp * D)
                ds_dbe = sig_unc / D

                do t = 1, flx_nobs
                    z = flx_obs(t) / sig
                    call dist_nll_score(flx_dist, z, p_dist1, p_dist2, nll_t, dl_dz, dp1, dp2)
                    factor = (1.0_dp + dl_dz * z) / sig
                    f            = f            + nll_t + log(sig)
                    grad_proc(1) = grad_proc(1) + factor * ds_dom
                    grad_proc(2) = grad_proc(2) + factor * ds_dal
                    grad_proc(3) = grad_proc(3) + factor * ds_dga
                    grad_proc(4) = grad_proc(4) + factor * ds_dbe
                    grad_dist1   = grad_dist1   + dp1
                    grad_dist2   = grad_dist2   + dp2

                    ! Advance recurrences: sigma_{t+1} = omega + (alpha+gamma*ind)*|y_t| + beta*sigma_t
                    if (flx_obs(t) < 0.0_dp) then
                        ind = 1.0_dp
                    else
                        ind = 0.0_dp
                    end if
                    ds_dom = 1.0_dp                         + beta * ds_dom
                    ds_dal = abs(flx_obs(t))                + beta * ds_dal
                    ds_dga = ind * abs(flx_obs(t))          + beta * ds_dga
                    ds_dbe = sig                            + beta * ds_dbe
                    sig    = omega + (alpha + gamma*ind) * abs(flx_obs(t)) + beta * sig
                end do

                ! Chain rule: same 3-way softmax Jacobian as GJR
                g(1) =  grad_proc(1) * omega
                g(2) =  alpha * ( grad_proc(2)*(1.0_dp - alpha) &
                                - 2.0_dp*grad_proc(3)*g2 &
                                - grad_proc(4)*beta )
                g(3) =  g2    * ( -grad_proc(2)*alpha &
                                + 2.0_dp*grad_proc(3)*(1.0_dp - g2) &
                                - grad_proc(4)*beta )
                g(4) =  beta  * ( -grad_proc(2)*alpha &
                                - 2.0_dp*grad_proc(3)*g2 &
                                + grad_proc(4)*(1.0_dp - beta) )
            end block

        ! ----------------------------------------------------------------
        case (proc_cnagarch)
            ! Two-component variance (Engle-Lee 1999):
            !   q_t = omega + rho*q_{t-1} + phi*(r²_{t-1} - h_{t-1})
            !   h_t = q_t + alpha*(r²_{t-1} - q_{t-1}) + beta*(h_{t-1} - q_{t-1})
            !   r_t = y_t - theta*sqrt(h_t),  z_t = y_t/sqrt(h_t)
            ! Structural params: omega(1), alpha(2), beta(3), rho(4), phi(5), theta(6)
            block
                real(dp) :: omega, alpha, beta, rho, phi, theta
                real(dp) :: h, q, sqrth, r, r2, C1, C2
                real(dp) :: A, B, kq
                real(dp) :: dh(6), dq(6), dh_new(6), dq_new(6)
                real(dp) :: dir_q(6), dir_h(6)
                integer  :: k
                real(dp), parameter :: eps_h = 1.0e-12_dp

                call cnagarch_transform(p(1:6), omega, alpha, beta, rho, phi, theta)

                ! Stationary initialisation: q_0 = omega/(1-rho), h_0 = q_0
                q = omega / max(1.0_dp - rho, eps_h)
                h = q

                ! Initialise derivatives of (h_0, q_0) w.r.t. structural params
                dq = 0.0_dp;  dh = 0.0_dp
                dq(1) = 1.0_dp / max(1.0_dp - rho, eps_h)   ! dq_0/domega
                dq(4) = q / max(1.0_dp - rho, eps_h)         ! dq_0/drho = omega/(1-rho)^2
                dh(1) = dq(1);  dh(4) = dq(4)

                do t = 1, flx_nobs
                    sqrth = sqrt(max(h, eps_h))
                    z = flx_obs(t) / sqrth
                    call dist_nll_score(flx_dist, z, p_dist1, p_dist2, nll_t, dl_dz, dp1, dp2)
                    factor = (1.0_dp + dl_dz * z) / (2.0_dp * max(h, eps_h))
                    f          = f          + nll_t + 0.5_dp * log(max(h, eps_h))
                    do k = 1, 6
                        grad_proc(k) = grad_proc(k) + factor * dh(k)
                    end do
                    grad_dist1 = grad_dist1 + dp1
                    grad_dist2 = grad_dist2 + dp2

                    ! Quantities for gradient recursion (use old h, q, sqrth)
                    r  = flx_obs(t) - theta * sqrth
                    r2 = r * r
                    C1 = r2 - q
                    C2 = h  - q
                    A  = rho - alpha - beta
                    B  = beta - (alpha + phi) * r * theta / sqrth - phi
                    kq = -phi * (1.0_dp + r * theta / sqrth)

                    ! Direct contributions per structural parameter
                    dir_q(1) = 1.0_dp;                     dir_h(1) = 0.0_dp
                    dir_q(2) = 0.0_dp;                     dir_h(2) = C1
                    dir_q(3) = 0.0_dp;                     dir_h(3) = C2
                    dir_q(4) = q;                          dir_h(4) = q
                    dir_q(5) = r2 - h;                     dir_h(5) = r2 - h
                    dir_q(6) = -2.0_dp*phi*r*sqrth;        dir_h(6) = -2.0_dp*(alpha+phi)*r*sqrth

                    do k = 1, 6
                        dq_new(k) = dir_q(k) + rho * dq(k) + kq  * dh(k)
                        dh_new(k) = dir_h(k) + A   * dq(k) + B   * dh(k)
                    end do
                    dq = dq_new
                    dh = dh_new

                    ! Advance state (C1, C2 used old q, h; update q first, then h)
                    q = omega + rho*q + phi*(r2 - h)
                    h = q + alpha*C1 + beta*C2
                    q = max(q, eps_h);  h = max(h, eps_h)
                end do

                ! Chain rule: unconstrained p -> structural params
                g(1) = grad_proc(1) * omega
                g(2) = grad_proc(2) * alpha*(1.0_dp - alpha) &
                     + grad_proc(3) * (-alpha * beta)
                g(3) = grad_proc(2) * (-alpha * beta) &
                     + grad_proc(3) * beta*(1.0_dp - beta)
                g(4) = grad_proc(4) * rho * (1.0_dp - rho)
                g(5) = grad_proc(5)
                g(6) = grad_proc(6)
            end block

        end select

        ! Distribution parameter gradients
        if (nd > 0) g(npp+1) = grad_dist1 / flx_nobs
        if (nd > 1) g(npp+2) = grad_dist2 / flx_nobs

        ! Normalise
        f         = f         / flx_nobs
        g(1:npp)  = g(1:npp)  / flx_nobs

    end subroutine flex_obj

    subroutine flex_skew_kurt(p, np, skew, kurt)
        ! Skew and excess kurtosis of standardised residuals z_t = y_t / sqrt(h_t).
        integer,  intent(in)  :: np
        real(dp), intent(in)  :: p(np)
        real(dp), intent(out) :: skew, kurt

        real(dp), allocatable :: zz(:)
        real(dp) :: zm, zv, zs, zk, dz, dz2, rn
        integer  :: t

        allocate(zz(flx_nobs))

        select case (flx_proc)
        case (proc_garch)
            block
                real(dp) :: omega, alpha, beta, h
                call garch_transform(p(1:3), omega, alpha, beta)
                h = omega / (1.0_dp - alpha - beta)
                do t = 1, flx_nobs
                    zz(t) = flx_obs(t) / sqrt(h)
                    h = omega + alpha*flx_obs(t)**2 + beta*h
                end do
            end block
        case (proc_nagarch)
            block
                real(dp) :: omega, alpha, beta, theta, h, s2, r
                call nagarch_transform(p(1:4), omega, alpha, beta, theta)
                call nagarch_shift_moments(theta, s2, r)
                h  = omega / (1.0_dp - alpha*s2 - beta)
                do t = 1, flx_nobs
                    zz(t) = flx_obs(t) / sqrt(h)
                    r = flx_obs(t) - theta*sqrt(h)
                    if (flx_nagarch_zero_above_shift .and. r > 0.0_dp) r = 0.0_dp
                    h = omega + alpha*r**2 + beta*h
                end do
            end block
        case (proc_gjr)
            block
                real(dp) :: omega, alpha, gamma, beta, h, g2, ind
                call gjr_transform(p(1:4), omega, alpha, gamma, beta)
                g2 = 0.5_dp * gamma
                h  = omega / (1.0_dp - alpha - g2 - beta)
                do t = 1, flx_nobs
                    zz(t) = flx_obs(t) / sqrt(h)
                    ind = merge(1.0_dp, 0.0_dp, flx_obs(t) < 0.0_dp)
                    h = omega + (alpha + gamma*ind)*flx_obs(t)**2 + beta*h
                end do
            end block
        case (proc_egarch)
            block
                real(dp) :: omega, alpha, gamma, beta, lh, z, c_eg
                call egarch_transform(p(1:4), omega, alpha, gamma, beta)
                c_eg = sqrt(2.0_dp / pi)
                lh   = omega / (1.0_dp - beta)
                do t = 1, flx_nobs
                    z     = flx_obs(t) / exp(0.5_dp*lh)
                    zz(t) = z
                    lh    = omega + beta*lh + alpha*(abs(z) - c_eg) + gamma*z
                end do
            end block
        case (proc_avgarch)
            block
                real(dp) :: omega, alpha, beta, sig
                call garch_transform(p(1:3), omega, alpha, beta)
                sig = omega / (1.0_dp - alpha - beta)
                do t = 1, flx_nobs
                    zz(t) = flx_obs(t) / sig
                    sig   = omega + alpha * abs(flx_obs(t)) + beta * sig
                end do
            end block
        case (proc_anagarch)
            block
                real(dp) :: omega, alpha, beta, theta, sig, r
                call anagarch_transform(p(1:4), omega, alpha, beta, theta)
                sig = omega / (1.0_dp - alpha - beta)
                do t = 1, flx_nobs
                    zz(t) = flx_obs(t) / sig
                    r     = flx_obs(t) - theta * sig
                    sig   = omega + alpha * abs(r) + beta * sig
                end do
            end block
        case (proc_agjr)
            block
                real(dp) :: omega, alpha, gamma, beta, sig, g2, ind
                call gjr_transform(p(1:4), omega, alpha, gamma, beta)
                g2  = 0.5_dp * gamma
                sig = omega / (1.0_dp - alpha - g2 - beta)
                do t = 1, flx_nobs
                    zz(t) = flx_obs(t) / sig
                    ind   = merge(1.0_dp, 0.0_dp, flx_obs(t) < 0.0_dp)
                    sig   = omega + (alpha + gamma*ind) * abs(flx_obs(t)) + beta * sig
                end do
            end block
        case (proc_cnagarch)
            block
                real(dp) :: omega, alpha, beta, rho, phi, theta
                real(dp) :: h, q, h_new, q_new, sqrth, r, r2
                call cnagarch_transform(p(1:6), omega, alpha, beta, rho, phi, theta)
                q = omega / max(1.0_dp - rho, 1.0e-12_dp)
                h = q
                do t = 1, flx_nobs
                    sqrth  = sqrt(max(h, 1.0e-12_dp))
                    zz(t)  = flx_obs(t) / sqrth
                    r      = flx_obs(t) - theta * sqrth
                    r2     = r * r
                    q_new  = omega + rho*q + phi*(r2 - h)
                    h_new  = q_new + alpha*(r2 - q) + beta*(h - q)
                    q = max(q_new, 1.0e-12_dp);  h = max(h_new, 1.0e-12_dp)
                end do
            end block
        end select

        rn = real(flx_nobs, dp)
        zm = sum(zz) / rn
        zv = 0.0_dp;  zs = 0.0_dp;  zk = 0.0_dp
        do t = 1, flx_nobs
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
    end subroutine flex_skew_kurt

    ! ----------------------------------------------------------------
    ! Component NAGARCH transform: p(6) -> (omega, alpha, beta, rho, phi, theta).
    ! omega = exp(p1); alpha, beta via 3-way softmax; rho = sigmoid(p4); phi, theta free.
    ! ----------------------------------------------------------------
    subroutine cnagarch_transform(p, omega, alpha, beta, rho, phi, theta)
        real(dp), intent(in)  :: p(6)
        real(dp), intent(out) :: omega, alpha, beta, rho, phi, theta
        real(dp) :: e2, e3, s
        omega = exp(p(1))
        e2    = exp(p(2));  e3 = exp(p(3))
        s     = 1.0_dp + e2 + e3
        alpha = e2 / s;  beta = e3 / s
        rho   = 1.0_dp / (1.0_dp + exp(-p(4)))
        phi   = p(5)
        theta = p(6)
    end subroutine cnagarch_transform

    subroutine cnagarch_inv_transform(omega, alpha, beta, rho, phi, theta, p)
        ! Requires omega > 0, alpha > 0, beta > 0, alpha+beta < 1, 0 < rho < 1.
        real(dp), intent(in)  :: omega, alpha, beta, rho, phi, theta
        real(dp), intent(out) :: p(6)
        real(dp) :: gam
        gam  = 1.0_dp - alpha - beta
        p(1) = log(omega)
        p(2) = log(alpha / gam)
        p(3) = log(beta  / gam)
        p(4) = log(rho / (1.0_dp - rho))
        p(5) = phi
        p(6) = theta
    end subroutine cnagarch_inv_transform

    ! ----------------------------------------------------------------
    ! Abs-NAGARCH transform: p(4) -> (omega, alpha, beta, theta).
    ! Parameterisation: same 3-param softmax as GARCH for (alpha, beta),
    ! so alpha+beta < 1 is enforced; theta (shift) is free (p4).
    ! ----------------------------------------------------------------
    subroutine anagarch_transform(p, omega, alpha, beta, theta)
        real(dp), intent(in)  :: p(4)
        real(dp), intent(out) :: omega, alpha, beta, theta  ! ANAGARCH parameters
        real(dp) :: e2, e3, s
        omega = exp(p(1))
        theta = p(4)
        e2    = exp(p(2))
        e3    = exp(p(3))
        s     = 1.0_dp + e2 + e3
        alpha = e2 / s
        beta  = e3 / s
    end subroutine anagarch_transform

    subroutine anagarch_inv_transform(omega, alpha, beta, theta, p)
        ! Requires omega > 0, alpha > 0, beta > 0, alpha+beta < 1.
        real(dp), intent(in)  :: omega, alpha, beta, theta  ! ANAGARCH parameters
        real(dp), intent(out) :: p(4)
        real(dp) :: gam
        gam  = 1.0_dp - alpha - beta
        p(1) = log(omega)
        p(2) = log(alpha / gam)
        p(3) = log(beta  / gam)
        p(4) = theta
    end subroutine anagarch_inv_transform

    ! ----------------------------------------------------------------
    ! Private: per-observation NLL from distribution and its scores.
    !
    !   nll_t  = -log f_std(z)   (distribution NLL, NO 0.5*log(h) term)
    !   dl_dz  = d(log f_std(z))/dz  (score, NOT negative)
    !   dp1    = d(nll_t)/d(p1)  (0 for distributions with no extra param)
    ! ----------------------------------------------------------------
    subroutine dist_nll_score(dist_type, z, p1, p2, nll_t, dl_dz, dp1, dp2)
        integer,  intent(in)  :: dist_type
        real(dp), intent(in)  :: z, p1, p2
        real(dp), intent(out) :: nll_t, dl_dz, dp1, dp2

        dp1 = 0.0_dp
        dp2 = 0.0_dp

        select case (dist_type)

        ! ---- Normal ----
        case (dist_normal)
            nll_t = log_sqrt_2pi + 0.5_dp * z**2
            dl_dz = -z

        ! ---- Student-t ----
        case (dist_t)
            block
                real(dp) :: nu, u, log1pu, f_const, dnu_const, dnll_dnu
                nu      = 2.0_dp + 98.0_dp / (1.0_dp + exp(-p1))
                u       = z**2 / (nu - 2.0_dp)
                log1pu  = log(1.0_dp + u)
                f_const = -log_gamma(0.5_dp*(nu+1.0_dp)) + log_gamma(0.5_dp*nu) &
                          + 0.5_dp * log(pi * (nu - 2.0_dp))
                nll_t   = f_const + 0.5_dp * (nu + 1.0_dp) * log1pu
                dl_dz   = -(nu + 1.0_dp) * z / (nu - 2.0_dp + z**2)
                dnu_const = -0.5_dp * digamma(0.5_dp*(nu+1.0_dp)) &
                            + 0.5_dp * digamma(0.5_dp*nu) &
                            + 0.5_dp / (nu - 2.0_dp)
                dnll_dnu = dnu_const + 0.5_dp*log1pu &
                           - 0.5_dp*(nu+1.0_dp)*u / ((nu-2.0_dp)*(1.0_dp+u))
                dp1 = dnll_dnu * (nu - 2.0_dp) * (100.0_dp - nu) / 98.0_dp
            end block

        ! ---- Hyperbolic Secant ----
        case (dist_sech)
            block
                real(dp) :: u
                u     = 0.5_dp * pi * z
                nll_t = log2 + log(cosh(u))
                dl_dz = -(0.5_dp * pi) * tanh(u)
            end block

        ! ---- GED ----
        case (dist_ged)
            block
                real(dp) :: nu, lam, a, a_nu, f_const, dfc_dnu, c_nu, dnll_dnu
                nu      = exp(p1)
                lam     = ged_lambda(nu)
                a       = abs(z) / lam
                if (a < tiny(1.0_dp)) a = tiny(1.0_dp)
                a_nu    = a**nu
                f_const = -log(nu) + log(2.0_dp) &
                          + 1.5_dp*log_gamma(1.0_dp/nu) - 0.5_dp*log_gamma(3.0_dp/nu)
                nll_t   = f_const + 0.5_dp * a_nu
                dl_dz   = -0.5_dp * nu * a_nu * sign(1.0_dp, z) &
                          / max(abs(z), 1.0e-300_dp)
                dfc_dnu = -1.0_dp/nu &
                          + (1.5_dp/nu**2) * (digamma(3.0_dp/nu) - digamma(1.0_dp/nu))
                c_nu    = (log(2.0_dp) - 0.5_dp*digamma(1.0_dp/nu) &
                          + 1.5_dp*digamma(3.0_dp/nu)) / nu
                if (a > 0.0_dp) then
                    dnll_dnu = dfc_dnu + 0.5_dp * a_nu * (log(a) - c_nu)
                else
                    dnll_dnu = dfc_dnu
                end if
                dp1 = dnll_dnu * nu
            end block

        ! ---- Laplace ----
        case (dist_laplace)
            nll_t = half_log2 + sqrt2 * abs(z)
            dl_dz = -sqrt2 * sign(1.0_dp, z)

        ! ---- Logistic ----
        case (dist_logistic)
            block
                real(dp) :: c_logis_loc, abs_u
                c_logis_loc = pi / sqrt3
                abs_u = abs(c_logis_loc * z)
                nll_t = -log(c_logis_loc) + abs_u + 2.0_dp * log(1.0_dp + exp(-abs_u))
                dl_dz = -c_logis_loc * sign(1.0_dp, z) * tanh(0.5_dp * abs_u)
            end block

        ! ---- NIG ----
        case (dist_nig)
            block
                real(dp), parameter :: alp_min = 0.1_dp
                real(dp), parameter :: alp_max = 20.0_dp
                real(dp), parameter :: alp_rng = 19.9_dp
                real(dp) :: alp, dalpd4, a2pu2, w, r, lk1, ratio
                real(dp) :: dnll_dalp
                alp    = alp_min + alp_rng / (1.0_dp + exp(-p1))
                dalpd4 = (alp - alp_min) * (alp_max - alp) / alp_rng
                a2pu2  = alp**2 + z**2
                w      = sqrt(a2pu2)
                r      = alp * w
                call bessel_k01(r, lk1, ratio)
                nll_t  = (log(pi) - 2.0_dp*log(alp) - alp**2) &
                         + 0.5_dp*log(a2pu2) - lk1
                dl_dz  = -z * (2.0_dp + r*ratio) / a2pu2
                dnll_dalp = (-1.0_dp/alp - 2.0_dp*alp) &
                            + ratio*(2.0_dp*alp**2 + z**2)/w &
                            + 2.0_dp*alp/a2pu2
                dp1 = dnll_dalp * dalpd4
            end block

        ! ---- Fernandez-Steel skewed-t ----
        case (dist_skew_t)
            block
                real(dp) :: nu, gamma, nll_p, nll_m, dl_dummy
                real(dp), parameter :: h_fd = 1.0e-5_dp

                nu    = 2.0_dp + 98.0_dp / (1.0_dp + exp(-p1))
                gamma = exp(p2)

                call skewt_nll_dl(z, nu, gamma, nll_t, dl_dz)

                ! FD gradient for p1 (p_nu, central differences)
                call skewt_nll_dl(z, 2.0_dp+98.0_dp/(1.0_dp+exp(-(p1+h_fd))), gamma, nll_p, dl_dummy)
                call skewt_nll_dl(z, 2.0_dp+98.0_dp/(1.0_dp+exp(-(p1-h_fd))), gamma, nll_m, dl_dummy)
                dp1 = (nll_p - nll_m) / (2.0_dp * h_fd)

                ! FD gradient for p2 (log gamma, central differences)
                call skewt_nll_dl(z, nu, gamma*exp( h_fd), nll_p, dl_dummy)
                call skewt_nll_dl(z, nu, gamma*exp(-h_fd), nll_m, dl_dummy)
                dp2 = (nll_p - nll_m) / (2.0_dp * h_fd)
            end block

        end select

    end subroutine dist_nll_score

    subroutine skewt_nll_dl(z, nu, gamma, nll, dl_dz)
        ! NLL and d(log f)/dz for the Fernandez-Steel skewed-t.
        ! Standardised so mean=0, variance=1 when nu>2.
        ! gamma>0 is the skew parameter (gamma=1 → symmetric t(nu)).
        real(dp), intent(in)  :: z, nu, gamma
        real(dp), intent(out) :: nll, dl_dz

        real(dp) :: c, mu1, sig2, sig, x, u, du_dz, log1pu2

        ! Peak of standardised t(nu): f_t(0) = Gamma((nu+1)/2) / (Gamma(nu/2)*sqrt(pi*(nu-2)))
        c = exp(log_gamma(0.5_dp*(nu+1.0_dp)) - log_gamma(0.5_dp*nu)) / sqrt(pi*(nu-2.0_dp))

        ! FS standardisation moments
        mu1  = 2.0_dp * (gamma - 1.0_dp/gamma) * c * (nu-2.0_dp) / (nu-1.0_dp)
        sig2 = gamma**2 + 1.0_dp/gamma**2 - 1.0_dp - mu1**2
        sig  = sqrt(sig2)

        ! Map to FS raw scale: x = z*sig + mu1
        x = z * sig + mu1

        if (x >= 0.0_dp) then
            u     = x / gamma
            du_dz = sig / gamma
        else
            u     = -x * gamma
            du_dz = -sig * gamma
        end if

        log1pu2 = log(1.0_dp + u**2 / (nu - 2.0_dp))

        ! NLL = -log f_FS(z; nu, gamma)
        nll = log(0.5_dp*(gamma + 1.0_dp/gamma)) &
              + log_gamma(0.5_dp*nu) - log_gamma(0.5_dp*(nu+1.0_dp)) &
              + 0.5_dp*log(pi*(nu-2.0_dp)) &
              + 0.5_dp*(nu+1.0_dp)*log1pu2 &
              + log(sig)

        ! d(log f)/dz = -(nu+1)*u/(nu-2+u^2) * du_dz
        dl_dz = -(nu+1.0_dp) * u / (nu-2.0_dp + u**2) * du_dz
    end subroutine skewt_nll_dl

    subroutine nagarch_shift_moments(theta, moment, dmoment)
        real(dp), intent(in)  :: theta
        real(dp), intent(out) :: moment, dmoment
        real(dp) :: cdf, pdf
        if (flx_nagarch_zero_above_shift) then
            cdf     = 0.5_dp * (1.0_dp + erf(theta / sqrt2))
            pdf     = exp(-0.5_dp * theta**2 - log_sqrt_2pi)
            moment  = (1.0_dp + theta**2) * cdf + theta * pdf
            dmoment = 2.0_dp * (theta * cdf + pdf)
        else
            moment  = 1.0_dp + theta**2
            dmoment = 2.0_dp * theta
        end if
    end subroutine nagarch_shift_moments

end module garch_flex_mod

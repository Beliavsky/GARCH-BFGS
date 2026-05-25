! Equal-weight ARCH: h_t = omega + alpha*(y^2_{t-1}+...+y^2_{t-k})/k
! Profiled over k = 1..K_max; GARCH(1,1)-Normal shown for comparison.
!
! Generalises historical volatility: omega=0, alpha=1 gives HV^2_k exactly.
! Only 2 free parameters (omega, alpha) regardless of k.
! NLL computed on observations k+1..n (neff = n-k); no initialisation needed.
!
! Parameters: p(1) = log(omega),  p(2) = logit(alpha)
!   omega > 0,  0 < alpha < 1  (stationarity: E[h] = omega/(1-alpha))

module arch_ew_mod
    use kind_mod,       only: dp
    use math_const_mod, only: log_sqrt_2pi
    implicit none
    private

    real(dp), allocatable, save :: ew_obs(:)
    real(dp), allocatable, save :: ew_csq(:)  ! cumulative sum of y_t^2, index 0..n
    integer,               save :: ew_nobs = 0
    integer,               save :: ew_k    = 1

    public :: ew_set_data, ew_set_k, ew_obj, ew_transform, ew_inv_transform

contains

    ! Store returns and precompute cumulative sum of squares.
    subroutine ew_set_data(y, n)
        integer,  intent(in) :: n     ! number of observations
        real(dp), intent(in) :: y(n)  ! demeaned returns
        integer :: t
        if (allocated(ew_obs)) deallocate(ew_obs)
        if (allocated(ew_csq)) deallocate(ew_csq)
        allocate(ew_obs(n), ew_csq(0:n))
        ew_obs    = y
        ew_nobs   = n
        ew_csq(0) = 0.0_dp
        do t = 1, n
            ew_csq(t) = ew_csq(t-1) + y(t)**2
        end do
    end subroutine ew_set_data

    ! Set window length before calling ew_obj.
    subroutine ew_set_k(k)
        integer, intent(in) :: k  ! lookback window (>= 1)
        ew_k = k
    end subroutine ew_set_k

    ! Unconstrained p(2) -> (omega, alpha).
    subroutine ew_transform(p, omega, alpha)
        real(dp), intent(in)  :: p(2)
        real(dp), intent(out) :: omega  ! level (> 0)
        real(dp), intent(out) :: alpha  ! weight on past k squared returns, in (0,1)
        omega = exp(p(1))
        alpha = 1.0_dp / (1.0_dp + exp(-p(2)))
    end subroutine ew_transform

    ! Constrained (omega, alpha) -> unconstrained p(2).
    subroutine ew_inv_transform(omega, alpha, p)
        real(dp), intent(in)  :: omega, alpha
        real(dp), intent(out) :: p(2)
        real(dp), parameter :: eps = 1.0e-10_dp
        real(dp) :: a
        p(1) = log(omega)
        a    = max(min(alpha, 1.0_dp - eps), eps)
        p(2) = log(a / (1.0_dp - a))
    end subroutine ew_inv_transform

    ! Normal NLL/neff and analytic gradient.
    ! s_t = (csq(t-1) - csq(t-1-k)) / k  computed in O(1) per step.
    subroutine ew_obj(p, np, f, g)
        integer,  intent(in)  :: np
        real(dp), intent(in)  :: p(np)
        real(dp), intent(out) :: f, g(np)
        real(dp) :: omega, alpha, h, s_t, fac
        real(dp) :: sum_domega, sum_dalpha, neff_r
        real(dp), parameter :: eps_h = 1.0e-12_dp
        integer :: t
        call ew_transform(p, omega, alpha)
        neff_r     = real(ew_nobs - ew_k, dp)
        f          = 0.0_dp
        sum_domega = 0.0_dp
        sum_dalpha = 0.0_dp
        do t = ew_k + 1, ew_nobs
            s_t = (ew_csq(t-1) - ew_csq(t-1-ew_k)) / real(ew_k, dp)
            h   = max(omega + alpha * s_t, eps_h)
            f   = f + log_sqrt_2pi + 0.5_dp*(log(h) + ew_obs(t)**2/h)
            fac = 0.5_dp * (1.0_dp - ew_obs(t)**2/h) / h
            sum_domega = sum_domega + fac
            sum_dalpha = sum_dalpha + fac * s_t
        end do
        f    = f          / neff_r
        g(1) = (sum_domega / neff_r) * omega                     ! d/d(log omega)
        g(2) = (sum_dalpha / neff_r) * alpha * (1.0_dp - alpha)  ! d/d(logit alpha)
    end subroutine ew_obj

end module arch_ew_mod


! ── Main program ─────────────────────────────────────────────────────────────
program xarch_ew
! Profile equal-weight ARCH over k = 1..K_max; compare against GARCH(1,1).

use kind_mod,       only: dp
use csv_mod,        only: read_price_csv, print_price_sample_info
use arch_ew_mod,    only: ew_set_data, ew_set_k, ew_obj, ew_transform, ew_inv_transform
use garch_mod,   only: garch_set_data, garch_obj, garch_transform, garch_inv_transform
use bfgs_mod,    only: bfgs_minimize
use stats_mod,      only: mean, sd
implicit none

integer,  parameter :: k_max        = 252     ! profile up to 1 trading year
real(dp), parameter :: trading_days = 252.0_dp
integer,  parameter :: max_iter     = 200
real(dp), parameter :: gtol         = 1.0e-6_dp
integer,  parameter :: nret         = 10**6

character(len=*), parameter :: prices_file = "vix_spy.csv"

integer,           allocatable :: dates(:)
character(len=32), allocatable :: col_names(:)
real(dp),          allocatable :: prices(:,:)
real(dp),          allocatable :: raw_ret(:), ret(:)
integer :: nobs, nprices, nall, ncols, icol, i1

real(dp) :: ret_mean, ret_std

real(dp), allocatable :: p(:), p0(:)
real(dp) :: fopt
integer  :: niter, k
logical  :: converged
real(dp) :: t_start, t_end

! EW-ARCH results
real(dp) :: ew_omega(k_max), ew_alpha(k_max), ew_vol(k_max)
real(dp) :: ew_logl_n(k_max), ew_aic(k_max), ew_bic(k_max)
integer  :: ew_neff(k_max)

! GARCH(1,1) results
real(dp) :: g_omega, g_alpha, g_beta, g_vol, g_logl_n, g_aic, g_bic

! working scalars
real(dp) :: omega, alpha, h_unc, vol_ann, logl, aic, bic, logl_n
integer  :: neff, k_aic, k_bic

! format strings: k(5), omega(12), alpha(8), vol%(9), neff(8), logL/n(9), AIC(12), BIC(12) = 75
character(len=*), parameter :: hdr_fmt = &
    "(A5,A12,A8,A9,A8,A9,A12,A12)"
character(len=*), parameter :: ew_row = &
    "(I5,ES12.3,F8.4,F9.2,I8,F9.4,F12.1,F12.1)"
character(len=*), parameter :: g_row  = &
    "(A5,ES12.3,F8.4,F9.2,I8,F9.4,F12.1,F12.1)"

call cpu_time(t_start)

call read_price_csv(prices_file, dates, col_names, prices)
nprices = size(prices, 1)
ncols   = size(col_names)
nall    = nprices - 1
nobs    = min(nret, nall)
i1      = nprices - nobs
allocate(raw_ret(nobs), ret(nobs))

call print_price_sample_info(prices_file, dates, ncols, nobs)
write(*, *)

do icol = 1, ncols

    raw_ret  = log(prices(i1+1:nprices,icol) / prices(i1:nprices-1,icol))
    ret_mean = mean(raw_ret)
    ret      = raw_ret - ret_mean
    ret_std  = sd(ret)

    write(*, '(A,A)') "Asset: ", trim(col_names(icol))
    write(*, '(A,F8.4,A,F8.4,A)') &
        "  mean (bps): ", ret_mean*1.0e4_dp, &
        "   emp vol_ann%: ", ret_std*sqrt(trading_days)*100.0_dp, "%"

    ! ── GARCH(1,1) baseline ──────────────────────────────────────────────────
    call garch_set_data(ret, nobs)
    allocate(p(3), p0(3))
    call garch_inv_transform(1.0e-5_dp, 0.08_dp, 0.88_dp, p0)
    p = p0
    call bfgs_minimize(garch_obj, p, 3, max_iter, gtol, fopt, niter, converged)
    call garch_transform(p, g_omega, g_alpha, g_beta)
    h_unc    = g_omega / (1.0_dp - g_alpha - g_beta)
    g_vol    = sqrt(trading_days * h_unc) * 100.0_dp
    g_logl_n = -fopt                                    ! log-likelihood per obs
    logl     = real(nobs,dp) * g_logl_n
    g_aic    = 2.0_dp*3  - 2.0_dp*logl
    g_bic    = 3*log(real(nobs,dp)) - 2.0_dp*logl
    deallocate(p, p0)

    ! ── EW-ARCH profile ───────────────────────────────────────────────────────
    call ew_set_data(ret, nobs)
    allocate(p(2), p0(2))
    call ew_inv_transform(1.0e-5_dp, 0.8_dp, p0)  ! common starting point for k=1

    do k = 1, k_max
        call ew_set_k(k)
        neff = nobs - k
        p = p0                          ! warm start from previous k's solution
        call bfgs_minimize(ew_obj, p, 2, max_iter, gtol, fopt, niter, converged)
        call ew_transform(p, omega, alpha)
        h_unc  = omega / max(1.0_dp - alpha, 1.0e-8_dp)
        vol_ann = sqrt(trading_days * h_unc) * 100.0_dp
        logl_n = -fopt
        logl   = real(neff,dp) * logl_n
        aic    = 4.0_dp - 2.0_dp*logl
        bic    = 2.0_dp*log(real(neff,dp)) - 2.0_dp*logl
        ew_omega(k)  = omega;   ew_alpha(k)  = alpha;  ew_vol(k)   = vol_ann
        ew_logl_n(k) = logl_n;  ew_aic(k)   = aic;    ew_bic(k)   = bic
        ew_neff(k)   = neff
        call ew_inv_transform(omega, alpha, p0)  ! warm start for k+1
    end do

    deallocate(p, p0)

    k_aic = minloc(ew_aic(1:k_max), 1)
    k_bic = minloc(ew_bic(1:k_max), 1)

    ! ── print ─────────────────────────────────────────────────────────────────
    write(*, hdr_fmt) &
        "    k", "       omega", "   alpha", " vol_ann%", "    neff", &
        "  logL/n", "         AIC", "         BIC"
    write(*, '(A)') repeat("-", 75)
    write(*, g_row) "GARCH", g_omega, g_alpha+g_beta, g_vol, nobs, g_logl_n, g_aic, g_bic
    write(*, '(A)') repeat("-", 75)

    do k = 1, k_max
        write(*, ew_row) k, ew_omega(k), ew_alpha(k), ew_vol(k), ew_neff(k), &
            ew_logl_n(k), ew_aic(k), ew_bic(k)
    end do

    write(*, '(A,I0,A,I0)') "  Best k: AIC -> ", k_aic, "   BIC -> ", k_bic
    write(*, '(A,F7.4,A,F7.4,A)') &
        "  logL/n: best EW-ARCH = ", ew_logl_n(k_aic), &
        "   GARCH(1,1) = ", g_logl_n, "  (directly comparable)"
    write(*, *)

end do

write(*, '(A)') "EW-ARCH: h_t = omega + alpha*(y^2_{t-1}+...+y^2_{t-k})/k   np=2"
write(*, '(A)') "GARCH(1,1): h_t = omega + alpha*y^2_{t-1} + beta*h_{t-1}   np=3"
write(*, '(A)') "alpha column: alpha for EW-ARCH; alpha+beta (persistence) for GARCH"
write(*, '(A)') "logL/n: per-observation log-likelihood — directly comparable across models"
write(*, '(A)') "AIC/BIC: within EW-ARCH for best-k selection (different n from GARCH)"

call cpu_time(t_end)
write(*, '(/,A,F0.3,A)') "elapsed time: ", t_end - t_start, " s"

end program xarch_ew

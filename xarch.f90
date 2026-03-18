! Unconstrained ARCH(k): h_t = omega + sum_{i=1}^{k} alpha_i * y^2_{t-i}
! All alpha_i > 0, sum(alpha_i) < 1 (stationarity); np = k+1 free parameters.
! Profiled over k = 1..K_max; GARCH(1,1)-Normal shown for comparison.
!
! Special case: k=1 is ARCH(1) with 2 parameters, identical to EW/LW at k=1.
! Unlike EW/LW, each lag has its own free coefficient -- np grows with k.
! NLL computed on observations k+1..n (neff = n-k); no initialisation needed.
!
! Parameters: p(1)   = log(omega)
!             p(i+1) = multinomial logit of alpha_i, i = 1..k
!             alpha_i = exp(p(i+1)) / (1 + sum_j exp(p(j+1)))
!             Ensures alpha_i > 0 and sum(alpha_i) < 1 for all k.
! Gradient via chain rule: g(j+1) = alpha_j * (G_j - sum_i alpha_i*G_i)
!   where G_i = dNLL/d(alpha_i) = (1/neff) * sum_t fac_t * y^2_{t-i}.

module arch_mod
    use kind_mod,       only: dp
    use math_const_mod, only: log_sqrt_2pi
    implicit none
    private

    real(dp), allocatable, save :: arch_obs(:)
    integer,               save :: arch_nobs = 0
    integer,               save :: arch_k    = 1

    public :: arch_set_data, arch_set_k, arch_obj, arch_transform, arch_inv_transform

contains

    ! Store returns.
    subroutine arch_set_data(y, n)
        integer,  intent(in) :: n     ! number of observations
        real(dp), intent(in) :: y(n)  ! demeaned returns
        if (allocated(arch_obs)) deallocate(arch_obs)
        allocate(arch_obs(n))
        arch_obs  = y
        arch_nobs = n
    end subroutine arch_set_data

    ! Set lag order before calling arch_obj.
    subroutine arch_set_k(k)
        integer, intent(in) :: k  ! lag order (>= 1)
        arch_k = k
    end subroutine arch_set_k

    ! Unconstrained p(np) -> omega and alpha(1:k), k = np-1.
    subroutine arch_transform(p, np, omega, alpha)
        integer,  intent(in)  :: np
        real(dp), intent(in)  :: p(np)
        real(dp), intent(out) :: omega
        real(dp), intent(out) :: alpha(np-1)   ! alpha(1:k), k = np-1
        real(dp) :: denom
        integer  :: i
        omega = exp(p(1))
        denom = 1.0_dp
        do i = 1, np-1
            denom = denom + exp(p(i+1))
        end do
        do i = 1, np-1
            alpha(i) = exp(p(i+1)) / denom
        end do
    end subroutine arch_transform

    ! Constrained (omega, alpha(1:k)) -> unconstrained p(k+1).
    subroutine arch_inv_transform(omega, alpha, k, p)
        real(dp), intent(in)  :: omega
        integer,  intent(in)  :: k
        real(dp), intent(in)  :: alpha(k)   ! lag coefficients
        real(dp), intent(out) :: p(k+1)
        real(dp), parameter :: eps = 1.0e-10_dp
        real(dp) :: leftover
        integer  :: i
        p(1)     = log(omega)
        leftover = max(1.0_dp - sum(alpha(1:k)), eps)
        do i = 1, k
            p(i+1) = log(max(alpha(i), eps) / leftover)
        end do
    end subroutine arch_inv_transform

    ! Normal NLL/neff and analytic gradient.
    ! np = k+1; p(1) = log(omega); p(2..np) = multinomial logit alphas.
    subroutine arch_obj(p, np, f, g)
        integer,  intent(in)  :: np
        real(dp), intent(in)  :: p(np)
        real(dp), intent(out) :: f, g(np)
        integer  :: k, t, i, j
        real(dp) :: omega, h, fac, neff_r, sum_ag
        real(dp) :: alpha(np-1), cap_g(np-1)
        real(dp), parameter :: eps_h = 1.0e-12_dp
        k = np - 1
        call arch_transform(p, np, omega, alpha)
        neff_r = real(arch_nobs - k, dp)
        f      = 0.0_dp
        g      = 0.0_dp
        cap_g  = 0.0_dp
        do t = k+1, arch_nobs
            h = omega
            do i = 1, k
                h = h + alpha(i) * arch_obs(t-i)**2
            end do
            h   = max(h, eps_h)
            f   = f + log_sqrt_2pi + 0.5_dp*(log(h) + arch_obs(t)**2/h)
            fac = 0.5_dp * (1.0_dp - arch_obs(t)**2/h) / h
            g(1) = g(1) + fac
            do i = 1, k
                cap_g(i) = cap_g(i) + fac * arch_obs(t-i)**2
            end do
        end do
        f     = f / neff_r
        g(1)  = (g(1) / neff_r) * omega           ! d/d(log omega)
        cap_g = cap_g / neff_r                     ! dNLL/d(alpha_i)
        sum_ag = dot_product(alpha, cap_g)
        do j = 1, k
            g(j+1) = alpha(j) * (cap_g(j) - sum_ag)   ! d/d(q_j)
        end do
    end subroutine arch_obj

end module arch_mod


! ── Main program ─────────────────────────────────────────────────────────────
program xarch
! Profile unconstrained ARCH(k) over k = 1..K_max; compare against GARCH(1,1).

use kind_mod,      only: dp
use csv_mod,       only: read_price_csv
use arch_mod,      only: arch_set_data, arch_set_k, arch_obj, arch_transform, arch_inv_transform
use garch_module,  only: garch_set_data, garch_obj, garch_transform, garch_inv_transform
use bfgs_module,   only: bfgs_minimize
use stats_mod,     only: mean, sd
implicit none

integer,  parameter :: k_max        = 20
real(dp), parameter :: trading_days = 252.0_dp
integer,  parameter :: max_iter     = 200
real(dp), parameter :: gtol         = 1.0e-6_dp
integer,  parameter :: nret         = 10**6

character(len=*), parameter :: prices_file = "vix_spy.csv"

integer,           allocatable :: dates(:)
character(len=32), allocatable :: col_names(:)
real(dp),          allocatable :: prices(:,:)
real(dp),          allocatable :: raw_ret(:), ret(:)
real(dp),          allocatable :: p_loc(:), p0_loc(:)
integer :: nobs, nprices, nall, ncols, icol, i1, i

real(dp) :: ret_mean, ret_std
real(dp) :: fopt
integer  :: niter, k, np
logical  :: converged
real(dp) :: t_start, t_end

! ARCH results
real(dp) :: ar_omega(k_max), ar_salpha(k_max), ar_vol(k_max)
real(dp) :: ar_logl_n(k_max), ar_aic(k_max), ar_bic(k_max)
integer  :: ar_neff(k_max)
real(dp) :: ar_alpha(k_max, k_max)   ! ar_alpha(k, 1:k) for ARCH(k)


! GARCH(1,1) results
real(dp) :: g_omega, g_alpha, g_beta, g_vol, g_logl_n, g_aic, g_bic

! working scalars
real(dp) :: omega, salpha, h_unc, vol_ann, logl, aic, bic, logl_n
real(dp) :: alpha_k(k_max), alpha_new(k_max), delta
real(dp) :: prev_omega, prev_alpha(k_max)
integer  :: neff, k_aic, k_bic

! format strings: k(5), omega(12), salpha(8), vol%(9), neff(8), logL/n(9), AIC(12), BIC(12) = 75
character(len=*), parameter :: hdr_fmt = &
    "(A5,A12,A8,A9,A8,A9,A12,A12)"
character(len=*), parameter :: ar_row = &
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

write(*, '(A,I0,A,I0,A)') "Using last ", nobs, " of ", nall, " observations"
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

    ! GARCH(1,1) baseline
    allocate(p_loc(3), p0_loc(3))
    call garch_set_data(ret, nobs)
    call garch_inv_transform(1.0e-5_dp, 0.08_dp, 0.88_dp, p0_loc)
    p_loc = p0_loc
    call bfgs_minimize(garch_obj, p_loc, 3, max_iter, gtol, fopt, niter, converged)
    call garch_transform(p_loc, g_omega, g_alpha, g_beta)
    h_unc    = g_omega / (1.0_dp - g_alpha - g_beta)
    g_vol    = sqrt(trading_days * h_unc) * 100.0_dp
    g_logl_n = -fopt
    logl     = real(nobs,dp) * g_logl_n
    g_aic    = 2.0_dp*3 - 2.0_dp*logl
    g_bic    = 3*log(real(nobs,dp)) - 2.0_dp*logl
    deallocate(p_loc, p0_loc)

    ! ARCH profile
    call arch_set_data(ret, nobs)
    ar_alpha   = 0.0_dp
    prev_alpha = 0.0_dp
    prev_omega = 0.0_dp

    do k = 1, k_max
        call arch_set_k(k)
        neff = nobs - k
        np   = k + 1
        allocate(p_loc(np), p0_loc(np))
        if (k == 1) then
            ! cold start: concentrate most weight in lag 1
            alpha_k(1) = 0.8_dp
            call arch_inv_transform(1.0e-5_dp, alpha_k(1:1), 1, p0_loc)
        else
            ! warm start: extend previous solution, trickle a small share to new lag
            delta = 0.05_dp
            alpha_new(1:k-1) = prev_alpha(1:k-1) * (1.0_dp - delta)
            alpha_new(k)     = delta * max(1.0_dp - sum(prev_alpha(1:k-1)), 1.0e-8_dp)
            call arch_inv_transform(prev_omega, alpha_new(1:k), k, p0_loc)
        end if
        p_loc = p0_loc
        call bfgs_minimize(arch_obj, p_loc, np, max_iter, gtol, fopt, niter, converged)
        call arch_transform(p_loc, np, omega, alpha_k(1:k))
        ar_omega(k)      = omega
        salpha           = sum(alpha_k(1:k))
        h_unc            = omega / max(1.0_dp - salpha, 1.0e-8_dp)
        vol_ann          = sqrt(trading_days * h_unc) * 100.0_dp
        logl_n           = -fopt
        logl             = real(neff,dp) * logl_n
        aic              = 2.0_dp*np - 2.0_dp*logl
        bic              = np*log(real(neff,dp)) - 2.0_dp*logl
        ar_salpha(k) = salpha;  ar_vol(k)    = vol_ann
        ar_logl_n(k) = logl_n; ar_aic(k)    = aic
        ar_bic(k)    = bic;    ar_neff(k)   = neff
        ar_alpha(k, 1:k) = alpha_k(1:k)
        prev_omega       = omega
        prev_alpha(1:k)  = alpha_k(1:k)
        deallocate(p_loc, p0_loc)
    end do

    k_aic = minloc(ar_aic(1:k_max), 1)
    k_bic = minloc(ar_bic(1:k_max), 1)

    write(*, hdr_fmt) &
        "    k", "       omega", "  salpha", " vol_ann%", "    neff", &
        "  logL/n", "         AIC", "         BIC"
    write(*, '(A)') repeat("-", 75)
    write(*, g_row) "GARCH", g_omega, g_alpha+g_beta, g_vol, nobs, g_logl_n, g_aic, g_bic
    write(*, '(A)') repeat("-", 75)

    do k = 1, k_max
        write(*, ar_row, advance='no') k, ar_omega(k), ar_salpha(k), ar_vol(k), ar_neff(k), &
            ar_logl_n(k), ar_aic(k), ar_bic(k)
        do i = 1, k
            write(*, '(F7.4)', advance='no') ar_alpha(k, i)
        end do
        write(*, *)
    end do

    write(*, '(A,I0,A,I0)') "  Best k: AIC -> ", k_aic, "   BIC -> ", k_bic
    write(*, '(A,F7.4,A,F7.4,A)') &
        "  logL/n: best ARCH = ", ar_logl_n(k_aic), &
        "   GARCH(1,1) = ", g_logl_n, "  (directly comparable)"
    write(*, *)

end do

write(*, '(A)') "ARCH(k): h_t = omega + sum_{i=1}^{k} alpha_i*y^2_{t-i}   np=k+1"
write(*, '(A)') "salpha = sum(alpha_i) -- total persistence"
write(*, '(A)') "GARCH(1,1): h_t = omega + alpha*y^2_{t-1} + beta*h_{t-1}   np=3"
write(*, '(A)') "salpha column: sum(alpha_i) for ARCH; alpha+beta (persistence) for GARCH"
write(*, '(A)') "logL/n: per-observation log-likelihood -- directly comparable across models"
write(*, '(A)') "AIC/BIC: within ARCH for best-k selection (different n from GARCH)"

call cpu_time(t_end)
write(*, '(/,A,F0.3,A)') "elapsed time: ", t_end - t_start, " s"

end program xarch

! ADCC-NAGARCH(1,1) Normal: two-stage estimation on N assets.
!
! Stage 1: fit NAGARCH(1,1)-Normal independently per asset.
!          h_t = omega + alpha*(r_{t-1} - theta*sqrt(h_{t-1}))^2 + beta*h_{t-1}
!          Unconditional variance: omega / (1 - alpha*(1+theta^2) - beta)
! Stage 2: fit ADCC(1,1) correlation on z_t = r_t/sqrt(h_t).
!          Q_t = (1-a-b)*Q_bar - g*N_bar + a*z_{t-1}z' + b*Q_{t-1} + g*n_{t-1}n'
!          where n_t = min(z_t, 0) captures asymmetric correlation responses.
!
! n_params = 4*N + 3.

program xadcc

use kind_mod,       only: dp
use csv_mod,        only: read_price_csv, print_price_sample_info
use nagarch_mod, only: nagarch_set_data, nagarch_obj, nagarch_transform, nagarch_inv_transform
use dcc_mod,        only: dcc_set_resid, adcc_obj, adcc_transform, adcc_inv_transform
use bfgs_mod,    only: bfgs_minimize
use stats_mod,      only: mean, sd
implicit none

real(dp), parameter :: trading_days = 252.0_dp
integer,  parameter :: max_iter     = 500
real(dp), parameter :: gtol         = 1.0e-6_dp
integer,  parameter :: nret         = 10**6

character(len=*), parameter :: prices_file = "spy_efa_eem_tlt_lqd.csv"

integer,           allocatable :: dates(:)
character(len=32), allocatable :: col_names(:)
real(dp),          allocatable :: prices(:,:)
real(dp),          allocatable :: ret_all(:,:)   ! (nobs, ncols)
real(dp),          allocatable :: h_all(:,:)     ! (nobs, ncols)
real(dp),          allocatable :: z_all(:,:)     ! (ncols, nobs)
real(dp),          allocatable :: p(:), p0(:)
integer :: nobs, nprices, nall, ncols, icol, i1

! Per-asset NAGARCH results
real(dp), allocatable :: g_omega(:), g_alpha(:), g_theta(:), g_beta(:)
real(dp), allocatable :: g_vol(:), g_logl_n(:)

! ADCC results
real(dp) :: adcc_a, adcc_b, adcc_g, fopt2
real(dp) :: logl_s1, logl_s2, total_logl, aic, bic

! working scalars
real(dp) :: omega, alpha_n, theta_n, beta_n, h_unc, fopt, ret_mean, persist
integer  :: niter, t, np_total
logical  :: converged
real(dp) :: t_start, t_end

! format strings
character(len=*), parameter :: hdr_fmt  = "(A10,A12,A8,A8,A8,A8,A9,A9)"
character(len=*), parameter :: asset_row = "(A10,ES12.3,F8.4,F8.4,F8.4,F8.4,F9.2,F9.4)"

call cpu_time(t_start)

call read_price_csv(prices_file, dates, col_names, prices)
nprices = size(prices, 1)
ncols   = size(col_names)
nall    = nprices - 1
nobs    = min(nret, nall)
i1      = nprices - nobs

allocate(ret_all(nobs,ncols), h_all(nobs,ncols), z_all(ncols,nobs))
allocate(g_omega(ncols), g_alpha(ncols), g_theta(ncols), g_beta(ncols))
allocate(g_vol(ncols), g_logl_n(ncols))

call print_price_sample_info(prices_file, dates, ncols, nobs)
write(*, *)

! ── Stage 1: NAGARCH(1,1) Normal per asset ───────────────────────────────────

allocate(p(4), p0(4))
do icol = 1, ncols
    ret_all(:,icol) = log(prices(i1+1:nprices,icol) / prices(i1:nprices-1,icol))
    ret_mean        = mean(ret_all(:,icol))
    ret_all(:,icol) = ret_all(:,icol) - ret_mean
    call nagarch_set_data(ret_all(:,icol), nobs)
    call nagarch_inv_transform(1.0e-5_dp, 0.08_dp, 0.88_dp, 0.0_dp, p0)
    p = p0
    call bfgs_minimize(nagarch_obj, p, 4, max_iter, gtol, fopt, niter, converged)
    call nagarch_transform(p, g_omega(icol), g_alpha(icol), g_beta(icol), g_theta(icol))
    g_logl_n(icol) = -fopt
    omega   = g_omega(icol);  alpha_n = g_alpha(icol)
    theta_n = g_theta(icol);  beta_n  = g_beta(icol)
    persist      = alpha_n*(1.0_dp + theta_n**2) + beta_n
    h_unc        = omega / max(1.0_dp - persist, 1.0e-8_dp)
    g_vol(icol)  = sqrt(trading_days * h_unc) * 100.0_dp
    ! Filter
    h_all(1,icol) = h_unc
    do t = 2, nobs
        h_all(t,icol) = omega + alpha_n*(ret_all(t-1,icol) - theta_n*sqrt(h_all(t-1,icol)))**2 &
                               + beta_n*h_all(t-1,icol)
    end do
    z_all(icol,:) = ret_all(:,icol) / sqrt(h_all(:,icol))
end do
deallocate(p, p0)

! Print stage-1 results
write(*, hdr_fmt) &
    "Asset     ", "       omega", "   alpha", "   theta", "    beta", " persist", " vol_ann%", "  logL1/n"
write(*, '(A)') repeat("-", 74)
do icol = 1, ncols
    persist = g_alpha(icol)*(1.0_dp + g_theta(icol)**2) + g_beta(icol)
    write(*, asset_row) col_names(icol), g_omega(icol), g_alpha(icol), g_theta(icol), &
        g_beta(icol), persist, g_vol(icol), g_logl_n(icol)
end do
write(*, *)

! ── Stage 2: ADCC(1,1) Normal ────────────────────────────────────────────────

call dcc_set_resid(z_all, ncols, nobs)
allocate(p(3), p0(3))
call adcc_inv_transform(0.05_dp, 0.88_dp, 0.05_dp, p0)
p = p0
call bfgs_minimize(adcc_obj, p, 3, max_iter, gtol, fopt2, niter, converged)
call adcc_transform(p, 3, adcc_a, adcc_b, adcc_g)
deallocate(p, p0)

! ── Totals ────────────────────────────────────────────────────────────────────

logl_s1    = sum(g_logl_n) * real(nobs, dp)
logl_s2    = -fopt2 * real(nobs, dp)
total_logl = logl_s1 + logl_s2
np_total   = ncols*4 + 3
aic        = 2.0_dp*np_total - 2.0_dp*total_logl
bic        = np_total*log(real(nobs,dp)) - 2.0_dp*total_logl

write(*, '(A,F8.4,A,F8.4,A,F8.4,A,F8.4)') &
    "ADCC: a = ", adcc_a, "   b = ", adcc_b, "   g = ", adcc_g, &
    "   a+b+g = ", adcc_a+adcc_b+adcc_g
write(*, '(A,I0)') "      converged = ", merge(1,0,converged)
write(*, *)
write(*, '(A,F10.4)') "Stage-1 logL/n   = ", logl_s1  / real(nobs, dp)
write(*, '(A,F10.4)') "Stage-2 logL/n   = ", logl_s2  / real(nobs, dp)
write(*, '(A,F10.4)') "Total   logL/n   = ", total_logl / real(nobs, dp)
write(*, '(A,I0)')    "n_params         = ", np_total
write(*, '(A,F12.1)') "AIC              = ", aic
write(*, '(A,F12.1)') "BIC              = ", bic
write(*, *)
write(*, '(A)') "Stage 1: NAGARCH(1,1)-Normal per asset, np=4 each"
write(*, '(A)') "Stage 2: ADCC(1,1)-Normal correlation, np=3"
write(*, '(A)') "persist = alpha*(1+theta^2) + beta"

call cpu_time(t_end)
write(*, '(/,A,F0.3,A)') "elapsed time: ", t_end - t_start, " s"

end program xadcc

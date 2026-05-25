! DCC-GARCH(1,1) Normal: two-stage estimation on N assets.
!
! Stage 1: fit GARCH(1,1)-Normal independently per asset.
! Stage 2: fit DCC(1,1) correlation on standardised residuals z_t = r_t/sqrt(h_t).
!
! Total log-likelihood = stage-1 logL (sum over assets) + stage-2 correlation logL.
! The stage-2 correction term [log|R_t| + z'R^{-1}z - z'z] exactly recovers the
! full joint Normal logL when combined with the marginal Normal stage-1 terms.
!
! AIC/BIC use n_params = 3*N + 2 (GARCH per asset + DCC a,b).

program xdcc

use kind_mod,      only: dp
use csv_mod,       only: read_price_csv, print_price_sample_info
use garch_mod,  only: garch_set_data, garch_obj, garch_transform, garch_inv_transform
use dcc_mod,       only: dcc_set_resid, dcc_obj, dcc_transform, dcc_inv_transform
use bfgs_mod,   only: bfgs_minimize
use stats_mod,     only: mean, sd
implicit none

real(dp), parameter :: trading_days = 252.0_dp
integer,  parameter :: max_iter     = 500
real(dp), parameter :: gtol         = 1.0e-6_dp
integer,  parameter :: nret         = 10**6

character(len=*), parameter :: prices_file = "spy_efa_eem_tlt_lqd.csv"

integer,           allocatable :: dates(:)
character(len=32), allocatable :: col_names(:)
real(dp),          allocatable :: prices(:,:)
real(dp),          allocatable :: ret_all(:,:)   ! (nobs, ncols) demeaned log returns
real(dp),          allocatable :: h_all(:,:)     ! (nobs, ncols) conditional variances
real(dp),          allocatable :: z_all(:,:)     ! (ncols, nobs) standardised residuals
real(dp),          allocatable :: p(:), p0(:)
integer :: nobs, nprices, nall, ncols, icol, i1

! Per-asset GARCH results
real(dp), allocatable :: g_omega(:), g_alpha(:), g_beta(:), g_vol(:), g_logl_n(:)

! DCC results
real(dp) :: dcc_a, dcc_b, fopt2
real(dp) :: logl_s1, logl_s2, total_logl, aic, bic

! working scalars
real(dp) :: omega, alpha_g, beta_g, h_unc, fopt, ret_mean
integer  :: niter, t, np_total
logical  :: converged
real(dp) :: t_start, t_end

! format strings
character(len=*), parameter :: hdr_fmt  = "(A10,A12,A8,A8,A8,A9,A9)"
character(len=*), parameter :: asset_row = "(A10,ES12.3,F8.4,F8.4,F8.4,F9.2,F9.4)"

call cpu_time(t_start)

call read_price_csv(prices_file, dates, col_names, prices)
nprices = size(prices, 1)
ncols   = size(col_names)
nall    = nprices - 1
nobs    = min(nret, nall)
i1      = nprices - nobs

allocate(ret_all(nobs,ncols), h_all(nobs,ncols), z_all(ncols,nobs))
allocate(g_omega(ncols), g_alpha(ncols), g_beta(ncols), g_vol(ncols), g_logl_n(ncols))

call print_price_sample_info(prices_file, dates, ncols, nobs)
write(*, *)

! ── Stage 1: GARCH(1,1) Normal per asset ─────────────────────────────────────

allocate(p(3), p0(3))
do icol = 1, ncols
    ret_all(:,icol) = log(prices(i1+1:nprices,icol) / prices(i1:nprices-1,icol))
    ret_mean        = mean(ret_all(:,icol))
    ret_all(:,icol) = ret_all(:,icol) - ret_mean
    call garch_set_data(ret_all(:,icol), nobs)
    call garch_inv_transform(1.0e-5_dp, 0.08_dp, 0.88_dp, p0)
    p = p0
    call bfgs_minimize(garch_obj, p, 3, max_iter, gtol, fopt, niter, converged)
    call garch_transform(p, g_omega(icol), g_alpha(icol), g_beta(icol))
    g_logl_n(icol) = -fopt
    alpha_g = g_alpha(icol);  beta_g = g_beta(icol);  omega = g_omega(icol)
    h_unc        = omega / max(1.0_dp - alpha_g - beta_g, 1.0e-8_dp)
    g_vol(icol)  = sqrt(trading_days * h_unc) * 100.0_dp
    ! Filter: h_t, z_t
    h_all(1,icol) = h_unc
    do t = 2, nobs
        h_all(t,icol) = omega + alpha_g*ret_all(t-1,icol)**2 + beta_g*h_all(t-1,icol)
    end do
    z_all(icol,:) = ret_all(:,icol) / sqrt(h_all(:,icol))
end do
deallocate(p, p0)

! Print stage-1 results
write(*, hdr_fmt) &
    "Asset     ", "       omega", "   alpha", "    beta", " persist", " vol_ann%", "  logL1/n"
write(*, '(A)') repeat("-", 66)
do icol = 1, ncols
    write(*, asset_row) col_names(icol), g_omega(icol), g_alpha(icol), g_beta(icol), &
        g_alpha(icol)+g_beta(icol), g_vol(icol), g_logl_n(icol)
end do
write(*, *)

! ── Stage 2: DCC(1,1) Normal ─────────────────────────────────────────────────

call dcc_set_resid(z_all, ncols, nobs)
allocate(p(2), p0(2))
call dcc_inv_transform(0.05_dp, 0.90_dp, p0)
p = p0
call bfgs_minimize(dcc_obj, p, 2, max_iter, gtol, fopt2, niter, converged)
call dcc_transform(p, 2, dcc_a, dcc_b)
deallocate(p, p0)

! ── Totals ────────────────────────────────────────────────────────────────────

logl_s1    = sum(g_logl_n) * real(nobs, dp)
logl_s2    = -fopt2 * real(nobs, dp)
total_logl = logl_s1 + logl_s2
np_total   = ncols*3 + 2
aic        = 2.0_dp*np_total - 2.0_dp*total_logl
bic        = np_total*log(real(nobs,dp)) - 2.0_dp*total_logl

write(*, '(A,F8.4,A,F8.4,A,F8.4)') &
    "DCC:  a = ", dcc_a, "   b = ", dcc_b, "   a+b = ", dcc_a+dcc_b
write(*, '(A,I0)') "      converged = ", merge(1,0,converged)
write(*, *)
write(*, '(A,F10.4)') "Stage-1 logL/n   = ", logl_s1  / real(nobs, dp)
write(*, '(A,F10.4)') "Stage-2 logL/n   = ", logl_s2  / real(nobs, dp)
write(*, '(A,F10.4)') "Total   logL/n   = ", total_logl / real(nobs, dp)
write(*, '(A,I0)')    "n_params         = ", np_total
write(*, '(A,F12.1)') "AIC              = ", aic
write(*, '(A,F12.1)') "BIC              = ", bic
write(*, *)
write(*, '(A)') "Stage 1: GARCH(1,1)-Normal per asset, np=3 each"
write(*, '(A)') "Stage 2: DCC(1,1)-Normal correlation, np=2"
write(*, '(A)') "logL/n: per-observation; stage-1 + stage-2 = full joint Normal logL/n"

call cpu_time(t_end)
write(*, '(/,A,F0.3,A)') "elapsed time: ", t_end - t_start, " s"

end program xdcc

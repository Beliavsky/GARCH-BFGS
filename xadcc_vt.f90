! ADCC-NAGARCH(1,1) Normal with variance targeting: two-stage estimation on N assets.
!
! Stage 1: NAGARCH(1,1)-VT per asset.  omega = sample_var*(1-alpha*(1+theta^2)-beta), np=3.
!          h_1 = sample_var exactly; unconditional variance anchored to sample.
! Stage 2: ADCC(1,1) Normal correlation.  Same as xadcc.f90.
!
! n_params = 3*N + 3.

program xadcc_vt

use kind_mod,       only: dp
use csv_mod,        only: read_price_csv
use nagarch_module, only: nagarch_set_data, nagarch_vt_set_target, nagarch_vt_obj, &
                          nagarch_vt_transform, nagarch_vt_inv_transform
use dcc_mod,        only: dcc_set_resid, adcc_obj, adcc_transform, adcc_inv_transform
use bfgs_module,    only: bfgs_minimize
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
real(dp),          allocatable :: ret_all(:,:)   ! (nobs, ncols) demeaned log returns
real(dp),          allocatable :: h_all(:,:)     ! (nobs, ncols) conditional variances
real(dp),          allocatable :: z_all(:,:)     ! (ncols, nobs) standardised residuals
real(dp),          allocatable :: p(:), p0(:)
integer :: nobs, nprices, nall, ncols, icol, i1

! Per-asset NAGARCH-VT results
real(dp), allocatable :: g_alpha(:), g_theta(:), g_beta(:), g_vol(:), g_logl_n(:)

! ADCC results
real(dp) :: adcc_a, adcc_b, adcc_g, fopt2
real(dp) :: logl_s1, logl_s2, total_logl, aic, bic

! working scalars
real(dp) :: alpha_n, theta_n, beta_n, omega, sample_var, persist, fopt, ret_mean
integer  :: niter, t, np_total
logical  :: converged
real(dp) :: t_start, t_end

! format strings
character(len=*), parameter :: hdr_fmt   = "(A10,A8,A8,A8,A8,A9,A9)"
character(len=*), parameter :: asset_row = "(A10,F8.4,F8.4,F8.4,F8.4,F9.2,F9.4)"

call cpu_time(t_start)

call read_price_csv(prices_file, dates, col_names, prices)
nprices = size(prices, 1)
ncols   = size(col_names)
nall    = nprices - 1
nobs    = min(nret, nall)
i1      = nprices - nobs

allocate(ret_all(nobs,ncols), h_all(nobs,ncols), z_all(ncols,nobs))
allocate(g_alpha(ncols), g_theta(ncols), g_beta(ncols), g_vol(ncols), g_logl_n(ncols))

write(*, '(A,I0,A,I0,A,I0,A)') &
    "Using last ", nobs, " of ", nall, " obs,  ", ncols, " assets"
write(*, *)

! ── Stage 1: NAGARCH(1,1)-VT per asset ───────────────────────────────────────

allocate(p(3), p0(3))
do icol = 1, ncols
    ret_all(:,icol) = log(prices(i1+1:nprices,icol) / prices(i1:nprices-1,icol))
    ret_mean        = mean(ret_all(:,icol))
    ret_all(:,icol) = ret_all(:,icol) - ret_mean
    sample_var      = sd(ret_all(:,icol))**2
    call nagarch_set_data(ret_all(:,icol), nobs)
    call nagarch_vt_set_target(sample_var)
    call nagarch_vt_inv_transform(0.08_dp, 0.88_dp, 0.0_dp, p0)
    p = p0
    call bfgs_minimize(nagarch_vt_obj, p, 3, max_iter, gtol, fopt, niter, converged)
    call nagarch_vt_transform(p, g_alpha(icol), g_beta(icol), g_theta(icol))
    g_logl_n(icol) = -fopt
    alpha_n = g_alpha(icol);  beta_n = g_beta(icol);  theta_n = g_theta(icol)
    omega         = sample_var * (1.0_dp - alpha_n*(1.0_dp + theta_n**2) - beta_n)
    persist       = alpha_n*(1.0_dp + theta_n**2) + beta_n
    g_vol(icol)   = sqrt(trading_days * sample_var) * 100.0_dp  ! = AnnVol% by construction
    ! Filter
    h_all(1,icol) = sample_var
    do t = 2, nobs
        h_all(t,icol) = omega + alpha_n*(ret_all(t-1,icol) - theta_n*sqrt(h_all(t-1,icol)))**2 &
                               + beta_n*h_all(t-1,icol)
    end do
    z_all(icol,:) = ret_all(:,icol) / sqrt(h_all(:,icol))
end do
deallocate(p, p0)

! Print stage-1 results
write(*, hdr_fmt) &
    "Asset     ", "   alpha", "   theta", "    beta", " persist", " vol_ann%", "  logL1/n"
write(*, '(A)') repeat("-", 66)
do icol = 1, ncols
    persist = g_alpha(icol)*(1.0_dp + g_theta(icol)**2) + g_beta(icol)
    write(*, asset_row) col_names(icol), g_alpha(icol), g_theta(icol), &
        g_beta(icol), persist, g_vol(icol), g_logl_n(icol)
end do
write(*, *)

! ── Stage 2: ADCC(1,1) Normal ─────────────────────────────────────────────────

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
np_total   = ncols*3 + 3
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
write(*, '(A)') "Stage 1: NAGARCH(1,1)-Normal-VT per asset (omega derived), np=3 each"
write(*, '(A)') "Stage 2: ADCC(1,1)-Normal correlation, np=3"
write(*, '(A)') "persist = alpha*(1+theta^2) + beta"
write(*, '(A)') "vol_ann% = sqrt(252*sample_var)*100 = AnnVol% by construction"

call cpu_time(t_end)
write(*, '(/,A,F0.3,A)') "elapsed time: ", t_end - t_start, " s"

end program xadcc_vt

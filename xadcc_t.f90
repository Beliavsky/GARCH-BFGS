! ADCC-NAGARCH(1,1) Student-t: two-stage estimation on N assets.
!
! Stage 1: fit NAGARCH(1,1)-Normal independently per asset (same as xadcc).
! Stage 2: fit ADCC(1,1) with multivariate Student-t noise.
!          Parameters: a, b, g (correlation dynamics) + nu > 2 (degrees of freedom).
!          NLL = (1/T)*sum_t [C_nu + 0.5*log|R_t| + 0.5*(nu+N)*log(1+z'R^{-1}z/(nu-2))]
!          where C_nu = log_gamma(nu/2) - log_gamma((nu+N)/2) + N/2*log(pi*(nu-2))
!
! The single shared nu ties all assets together: joint tail events are more likely
! than under independent marginals.  Compare AIC/BIC with xadcc to assess fit gain.
!
! n_params = 4*N + 4.

program xadcc_t

use kind_mod,       only: dp
use csv_mod,        only: read_price_csv, print_price_sample_info
use nagarch_mod, only: nagarch_set_data, nagarch_obj, nagarch_transform, nagarch_inv_transform
use dcc_mod,        only: dcc_set_resid, adcc_t_obj, adcc_t_transform, adcc_t_inv_transform
use bfgs_mod,    only: bfgs_minimize
use stats_mod,      only: mean, sd
implicit none

real(dp), parameter :: trading_days = 252.0_dp
integer,  parameter :: max_iter     = 500
real(dp), parameter :: gtol         = 1.0e-6_dp
integer,  parameter :: nret         = 10**6
integer,  parameter :: n_ac         = 5    ! autocorrelation lags for squared returns

character(len=*), parameter :: prices_file = "spy_efa_eem_tlt_lqd.csv"

integer,           allocatable :: dates(:)
character(len=32), allocatable :: col_names(:)
real(dp),          allocatable :: prices(:,:)
real(dp),          allocatable :: ret_all(:,:)   ! (nobs, ncols) demeaned log returns
real(dp),          allocatable :: h_all(:,:)     ! (nobs, ncols) conditional variances
real(dp),          allocatable :: z_all(:,:)     ! (ncols, nobs) standardised residuals
real(dp),          allocatable :: ret_means(:)   ! raw daily mean per asset (for ann. return)
real(dp),          allocatable :: sq_ret(:)      ! squared returns workspace
real(dp),          allocatable :: p(:), p0(:)
integer :: nobs, nprices, nall, ncols, icol, i1

! Per-asset NAGARCH results
real(dp), allocatable :: g_omega(:), g_alpha(:), g_theta(:), g_beta(:)
real(dp), allocatable :: g_vol(:), g_logl_n(:)

! ADCC-t results
real(dp) :: adcc_a, adcc_b, adcc_g, adcc_nu, fopt2
real(dp) :: logl_s1, logl_s2, total_logl, aic, bic, h_logsum

! working scalars
real(dp) :: omega, alpha_n, theta_n, beta_n, h_unc, fopt, persist
real(dp) :: ret_std, ann_ret, ann_vol, skew, ekurt, sq_mean, sq_var, ac_j
integer  :: niter, t, j, icol2, np_total
logical  :: converged
real(dp) :: t_start, t_end

! format strings
character(len=*), parameter :: hdr_fmt   = "(A10,A12,A8,A8,A8,A8,A9,A9)"
character(len=*), parameter :: asset_row = "(A10,ES12.3,F8.4,F8.4,F8.4,F8.4,F9.2,F9.4)"
character(len=*), parameter :: stats_hdr = "(A10,A9,A9,A9,A9)"
character(len=*), parameter :: stats_row = "(A10,F9.2,F9.2,F9.3,F9.3)"
character(len=*), parameter :: ac_hdr    = "(A9)"
character(len=*), parameter :: ac_val    = "(F9.4)"

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
allocate(ret_means(ncols), sq_ret(nobs))

call print_price_sample_info(prices_file, dates, ncols, nobs)
write(*, *)

! ── Compute log returns ───────────────────────────────────────────────────────

do icol = 1, ncols
    ret_all(:,icol) = log(prices(i1+1:nprices,icol) / prices(i1:nprices-1,icol))
    ret_means(icol) = mean(ret_all(:,icol))
    ret_all(:,icol) = ret_all(:,icol) - ret_means(icol)
end do

! ── Descriptive statistics ────────────────────────────────────────────────────

write(*, stats_hdr, advance='no') &
    "Asset     ", " AnnRet%", " AnnVol%", "    Skew", " ExKurt"
do j = 1, n_ac
    write(*, ac_hdr, advance='no') "      AC" // char(48+j)
end do
write(*, *)
write(*, '(A)') repeat("-", 10 + 9*4 + 9*n_ac)

do icol = 1, ncols
    ret_std = sd(ret_all(:,icol))
    ann_ret = ret_means(icol) * trading_days * 100.0_dp
    ann_vol = ret_std * sqrt(trading_days) * 100.0_dp

    ! Skewness and excess kurtosis (standardised moments)
    skew  = sum((ret_all(:,icol)/ret_std)**3) / real(nobs, dp)
    ekurt = sum((ret_all(:,icol)/ret_std)**4) / real(nobs, dp) - 3.0_dp

    ! Autocorrelations of squared returns
    sq_ret  = ret_all(:,icol)**2
    sq_mean = sum(sq_ret) / real(nobs, dp)
    sq_var  = sum((sq_ret - sq_mean)**2)   ! denominator (same for all lags)

    write(*, stats_row, advance='no') col_names(icol), ann_ret, ann_vol, skew, ekurt
    do j = 1, n_ac
        ac_j = sum((sq_ret(j+1:nobs) - sq_mean) * (sq_ret(1:nobs-j) - sq_mean)) / sq_var
        write(*, ac_val, advance='no') ac_j
    end do
    write(*, *)
end do
write(*, *)

! ── Stage 1: NAGARCH(1,1) Normal per asset ───────────────────────────────────

allocate(p(4), p0(4))
do icol = 1, ncols
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

! ── Stage 2: ADCC(1,1) Student-t ─────────────────────────────────────────────

call dcc_set_resid(z_all, ncols, nobs)
allocate(p(4), p0(4))
call adcc_t_inv_transform(0.05_dp, 0.88_dp, 0.05_dp, 8.0_dp, p0)
p = p0
call bfgs_minimize(adcc_t_obj, p, 4, max_iter, gtol, fopt2, niter, converged)
call adcc_t_transform(p, 4, adcc_a, adcc_b, adcc_g, adcc_nu)
deallocate(p, p0)

! ── Totals ────────────────────────────────────────────────────────────────────
! Full joint multivariate-t logL = stage-2 t logL - 0.5*sum_{i,t} log(h_it).
! The stage-2 t NLL already contains the log-Gamma and log|R_t| terms but not
! the marginal variance terms -0.5*log(h_it); adding them gives the full joint
! t log-likelihood, directly comparable to the Normal joint logL in xdcc/xadcc.

logl_s1  = sum(g_logl_n) * real(nobs, dp)   ! Normal stage-1 (for reference only)
logl_s2  = -fopt2 * real(nobs, dp)          ! t stage-2 (log-Gamma terms included)
h_logsum = 0.0_dp
do icol2 = 1, ncols
    do t = 1, nobs
        h_logsum = h_logsum + log(h_all(t, icol2))
    end do
end do
total_logl = logl_s2 - 0.5_dp * h_logsum    ! full joint t logL
np_total   = ncols*4 + 4
aic        = 2.0_dp*np_total - 2.0_dp*total_logl
bic        = np_total*log(real(nobs,dp)) - 2.0_dp*total_logl

write(*, '(A,F8.4,A,F8.4,A,F8.4,A,F7.3)') &
    "ADCC-t: a = ", adcc_a, "   b = ", adcc_b, "   g = ", adcc_g, &
    "   nu = ", adcc_nu
write(*, '(A,F8.4,A,I0)') &
    "        a+b+g = ", adcc_a+adcc_b+adcc_g, "   converged = ", merge(1,0,converged)
write(*, *)
write(*, '(A,F10.4)') "Stage-1 Normal logL/n = ", logl_s1  / real(nobs, dp)
write(*, '(A,F10.4)') "Stage-2 t      logL/n = ", logl_s2  / real(nobs, dp)
write(*, '(A,F10.4)') "Full joint t   logL/n = ", total_logl / real(nobs, dp)
write(*, '(A,I0)')    "n_params              = ", np_total
write(*, '(A,F12.1)') "AIC                   = ", aic
write(*, '(A,F12.1)') "BIC                   = ", bic
write(*, *)
write(*, '(A)') "Stage 1: NAGARCH(1,1)-Normal per asset, np=4 each"
write(*, '(A)') "Stage 2: ADCC(1,1)-Student-t correlation, np=4 (a, b, g, nu)"
write(*, '(A)') "persist = alpha*(1+theta^2) + beta"
write(*, '(A)') "Full joint t logL = stage-2 t logL - 0.5*sum log(h_it); comparable to Normal total logL"

call cpu_time(t_end)
write(*, '(/,A,F0.3,A)') "elapsed time: ", t_end - t_start, " s"

end program xadcc_t

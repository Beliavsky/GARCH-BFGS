! Fit supported standardised distributions to iid returns from a price CSV file.
!
! Stage 1: fit NAGARCH(1,1)-Normal per asset to extract standardised residuals z_t.
! Stage 2: fit supported distributions to each asset's z_t via fit_dist_std,
!          then to the raw returns via fit_dist.
!
! Output table: for each asset x distribution: logL/n and shape parameters.

program xfit_dist_returns
    use date_mod, only: print_program_header

use kind_mod,          only: dp
use csv_mod,           only: read_price_csv, print_price_sample_info
use nagarch_mod,    only: nagarch_set_data, nagarch_obj, nagarch_transform, nagarch_inv_transform
use bfgs_mod,       only: bfgs_minimize
use stats_mod,         only: mean, sd
use distributions_mod, only: dist_count, dist_names, dist_npar_std, &
                              fit_dist_std, fit_dist
implicit none

real(dp), parameter :: trading_days = 252.0_dp
integer,  parameter :: max_iter     = 500
real(dp), parameter :: gtol         = 1.0e-6_dp
integer,  parameter :: nret         = 10**6

character(len=*), parameter :: prices_file = "spy_efa_eem_tlt_lqd.csv"

integer,           allocatable :: dates(:)
character(len=32), allocatable :: col_names(:)
real(dp),          allocatable :: prices(:,:)
real(dp),          allocatable :: ret_all(:,:)
real(dp),          allocatable :: h_all(:,:)
real(dp),          allocatable :: z_all(:,:)
real(dp),          allocatable :: p(:), p0(:)
integer :: nobs, nprices, nall, ncols, icol, i1, d

real(dp) :: shape, loglik, mu, sigma
real(dp) :: omega, alpha_n, beta_n, theta_n, h_unc, persist, fopt, ret_mean
integer  :: niter, t
logical  :: converged, conv_d
real(dp) :: t_start, t_end

    call print_program_header("xfit_dist_returns.f90")
call cpu_time(t_start)

call read_price_csv(prices_file, dates, col_names, prices)
nprices = size(prices, 1)
ncols   = size(col_names)
nall    = nprices - 1
nobs    = min(nret, nall)
i1      = nprices - nobs

allocate(ret_all(nobs,ncols), h_all(nobs,ncols), z_all(ncols,nobs))

call print_price_sample_info(prices_file, dates, ncols, nobs)
write(*, *)

! ?????? Stage 1: NAGARCH(1,1)-Normal per asset ?????????????????????????????????????????????????????????????????????????????????????????????????????????

allocate(p(4), p0(4))
do icol = 1, ncols
    ret_all(:,icol) = log(prices(i1+1:nprices,icol) / prices(i1:nprices-1,icol))
    ret_mean        = mean(ret_all(:,icol))
    ret_all(:,icol) = ret_all(:,icol) - ret_mean
    call nagarch_set_data(ret_all(:,icol), nobs)
    call nagarch_inv_transform(1.0e-5_dp, 0.08_dp, 0.88_dp, 0.0_dp, p0)
    p = p0
    call bfgs_minimize(nagarch_obj, p, 4, max_iter, gtol, fopt, niter, converged)
    call nagarch_transform(p, omega, alpha_n, beta_n, theta_n)
    persist      = alpha_n*(1.0_dp + theta_n**2) + beta_n
    h_unc        = omega / max(1.0_dp - persist, 1.0e-8_dp)
    h_all(1,icol) = h_unc
    do t = 2, nobs
        h_all(t,icol) = omega + alpha_n*(ret_all(t-1,icol) - theta_n*sqrt(h_all(t-1,icol)))**2 &
                               + beta_n*h_all(t-1,icol)
    end do
    z_all(icol,:) = ret_all(:,icol) / sqrt(h_all(:,icol))
end do
deallocate(p, p0)

! ?????? Stage 2: fit distributions to standardised residuals z_t ???????????????????????????????????????????????????

write(*, '(A)') "Distribution fit to NAGARCH(1,1) standardised residuals z_t:"
write(*, '(A)') "logL/n and shape parameter for each asset x distribution"
write(*, *)

! Header
write(*, '(A10)', advance='no') "Dist      "
do icol = 1, ncols
    write(*, '(A14,A5)', advance='no') adjustr(col_names(icol)(1:8)), " shp "
end do
write(*, *)
write(*, '(A)') repeat("-", 10 + ncols*19)

do d = 1, dist_count
    write(*, '(A10)', advance='no') trim(dist_names(d))
    do icol = 1, ncols
        call fit_dist_std(z_all(icol,:), nobs, d, shape, loglik, conv_d)
        if (dist_npar_std(d) == 0) then
            write(*, '(F9.4,5X,A5)', advance='no') loglik/real(nobs,dp), "  -- "
        else
            write(*, '(F9.4,F9.3)', advance='no') loglik/real(nobs,dp), shape
        end if
    end do
    write(*, *)
end do

write(*, *)
write(*, '(A)') "shape: nu (t, ged), alp (nig); -- = no shape parameter"
write(*, '(A)') "npar_std: normal/logistic/laplace/sech=0, t/ged/nig=1"

! ?????? Also: fit distributions to raw demeaned returns (full location-scale) ????????????

write(*, *)
write(*, '(A)') "Distribution fit to raw demeaned returns (mu=sample mean, sigma+shape fitted):"
write(*, '(A)') "logL/n, sigma (annualised vol%), shape"
write(*, *)

write(*, '(A10)', advance='no') "Dist      "
do icol = 1, ncols
    write(*, '(A12,A7)', advance='no') adjustr(col_names(icol)(1:8)), " shp   "
end do
write(*, *)
write(*, '(A)') repeat("-", 10 + ncols*19)

do d = 1, dist_count
    write(*, '(A10)', advance='no') trim(dist_names(d))
    do icol = 1, ncols
        call fit_dist(ret_all(:,icol), nobs, d, mu, sigma, shape, loglik, conv_d)
        if (dist_npar_std(d) == 0) then
            write(*, '(F9.4,F5.1,A5)', advance='no') &
                loglik/real(nobs,dp), sigma*sqrt(trading_days)*100.0_dp, "  -- "
        else
            write(*, '(F9.4,F5.1,F4.1)', advance='no') &
                loglik/real(nobs,dp), sigma*sqrt(trading_days)*100.0_dp, shape
        end if
    end do
    write(*, *)
end do

write(*, *)
write(*, '(A)') "sigma column: annualised vol% = sigma*sqrt(252)*100"
write(*, '(A)') "shape: nu (t, ged), alp (nig); -- = no shape parameter"

call cpu_time(t_end)
write(*, '(/,A,F0.3,A)') "elapsed time: ", t_end - t_start, " s"

end program xfit_dist_returns

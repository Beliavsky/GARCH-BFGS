program xfit_garch_returns
! Read price data from a CSV and fit GARCH(1,1), NAGARCH(1,1),
! GJR-GARCH(1,1), and EGARCH(1,1) with Student-t innovations to the
! log-returns of each column.  Reports parameter estimates, annualised
! unconditional vol, log-likelihood, AIC, BIC, and model ranks.

use kind_mod,         only: dp
use csv_mod,          only: read_price_csv, print_price_sample_info
use garch_flex_mod,   only: flex_set_data, flex_set_types, flex_np, flex_obj, &
                             proc_garch, proc_nagarch, proc_gjr, proc_egarch, &
                             proc_names, dist_t
use garch_mod,     only: garch_inv_transform,   garch_transform
use nagarch_mod,   only: nagarch_inv_transform, nagarch_transform
use gjr_mod,       only: gjr_inv_transform,     gjr_transform
use egarch_mod,    only: egarch_inv_transform,  egarch_transform
use bfgs_mod,      only: bfgs_minimize
use stats_mod,        only: mean, sd
use rank_mod,         only: rank_desc, rank_asc
implicit none

integer, parameter :: nmod = proc_egarch   ! number of process types

! constants
real(dp), parameter :: trading_days = 252.0_dp
real(dp), parameter :: p_t0         = -2.729_dp  ! t-dist start: nu ≈ 8
integer,  parameter :: max_iter     = 100
real(dp), parameter :: gtol         = 1.0e-7_dp
integer,  parameter :: nret         = 10**6       ! cap on observations used

! data
integer,           allocatable :: dates(:)
character(len=32), allocatable :: col_names(:)
real(dp),          allocatable :: prices(:,:)
real(dp),          allocatable :: raw_ret(:)        ! raw log-returns before demeaning
real(dp),          allocatable :: ret(:)           ! demeaned log-returns
integer :: nobs, nprices, nall, ncols, icol

real(dp) :: ret_mean, ret_std

! optimisation
real(dp), allocatable :: p(:), p0(:)
real(dp) :: fopt
integer  :: np, niter, iproc
logical  :: converged
real(dp) :: t_start, t_end

! per-model results
real(dp) :: omegas(nmod), alphas(nmod), gamma_ps(nmod), betas(nmod)
real(dp) :: nus(nmod), vol_anns(nmod), logls(nmod), aics(nmod), bics(nmod)
integer  :: niters(nmod)
logical  :: conv(nmod)

! working scalars
real(dp) :: omega, alpha, gamma_p, beta, theta, nu
real(dp) :: h_unc, vol_ann, logl, aic, bic
character (len=*), parameter :: prices_file = "vix_spy.csv"
integer          :: rank_logl(nmod), rank_aic(nmod), rank_bic(nmod)

! read data once
call read_price_csv(prices_file, dates, col_names, prices)
nprices = size(prices, 1)
ncols   = size(col_names)
nall    = nprices - 1
nobs    = min(nret, nall)
allocate(raw_ret(nobs), ret(nobs))

call print_price_sample_info(prices_file, dates, ncols, nobs)
write(*,*)

call cpu_time(t_start)

! loop over assets
do icol = 1, ncols

    raw_ret  = log(prices(nall-nobs+2:nprices, icol) / prices(nall-nobs+1:nprices-1, icol))
    ret_mean = mean(raw_ret)
    ret      = raw_ret - ret_mean
    ret_std  = sd(ret)

    write(*, '(A,A)')    "Asset: ", trim(col_names(icol))
    write(*, '(A,F8.4,A,F8.4,A)') &
        "  mean (bps): ", ret_mean*1.0e4_dp, "   emp vol_ann%: ", ret_std*sqrt(trading_days)*100.0_dp, "%"
    write(*, '(A8,A12,A8,A9,A9,A7,A9,A12,A12,A12,A6,A8,A6,A6)') &
        "Process", "omega", "alpha", "par3", "beta", "nu", &
        "vol_ann%", "logL", "AIC", "BIC", "iter", "rklogL", "rkAIC", "rkBIC"
    write(*, '(A)') repeat("-", 127)

    do iproc = proc_garch, proc_egarch

        call flex_set_data(ret, nobs)
        call flex_set_types(iproc, dist_t)
        np = flex_np()
        allocate(p(np), p0(np))

        select case (iproc)
        case (proc_garch)
            call garch_inv_transform(1.0e-5_dp, 0.08_dp, 0.88_dp, p0(1:3))

        case (proc_nagarch)
            ! theta > 0: negative shocks amplify variance for equity indices
            call nagarch_inv_transform(1.0e-5_dp, 0.06_dp, 0.88_dp, 0.5_dp, p0(1:4))

        case (proc_gjr)
            call gjr_inv_transform(1.0e-5_dp, 0.04_dp, 0.08_dp, 0.88_dp, p0(1:4))

        case (proc_egarch)
            ! omega start reflects unconditional log-variance: log(ret_std^2)*(1-beta)
            call egarch_inv_transform(log(ret_std**2)*0.03_dp, 0.10_dp, -0.10_dp, 0.97_dp, p0(1:4))
        end select

        p0(np) = p_t0
        p = p0

        call bfgs_minimize(flex_obj, p, np, max_iter, gtol, fopt, niter, converged)

        nu = 2.0_dp + 98.0_dp / (1.0_dp + exp(-p(np)))

        select case (iproc)
        case (proc_garch)
            call garch_transform(p(1:3), omega, alpha, beta)
            gamma_p = 0.0_dp
            h_unc   = omega / (1.0_dp - alpha - beta)

        case (proc_nagarch)
            call nagarch_transform(p(1:4), omega, alpha, beta, theta)
            gamma_p = theta
            h_unc   = omega / (1.0_dp - alpha*(1.0_dp + theta**2) - beta)

        case (proc_gjr)
            call gjr_transform(p(1:4), omega, alpha, gamma_p, beta)
            h_unc   = omega / (1.0_dp - alpha - 0.5_dp*gamma_p - beta)

        case (proc_egarch)
            call egarch_transform(p(1:4), omega, alpha, gamma_p, beta)
            h_unc   = exp(omega / (1.0_dp - beta))
        end select

        vol_ann = sqrt(trading_days * h_unc) * 100.0_dp
        logl    = -nobs * fopt
        aic     = 2.0_dp*np             - 2.0_dp*logl
        bic     = np*log(real(nobs,dp)) - 2.0_dp*logl

        omegas(iproc)   = omega
        alphas(iproc)   = alpha
        gamma_ps(iproc) = gamma_p
        betas(iproc)    = beta
        nus(iproc)      = nu
        vol_anns(iproc) = vol_ann
        logls(iproc)    = logl
        aics(iproc)     = aic
        bics(iproc)     = bic
        niters(iproc)   = niter
        conv(iproc)     = converged

        deallocate(p, p0)

    end do

    call rank_desc(logls, rank_logl)
    call rank_asc(aics,   rank_aic)
    call rank_asc(bics,   rank_bic)

    do iproc = proc_garch, proc_egarch
        if (iproc == proc_garch) then
            write(*, '(A8,ES12.3,F8.4,9X,F9.4,F7.2,F9.2,F12.1,F12.1,F12.1,I6,3I6)') &
                trim(proc_names(iproc)), omegas(iproc), alphas(iproc), &
                betas(iproc), nus(iproc), vol_anns(iproc), &
                logls(iproc), aics(iproc), bics(iproc), niters(iproc), &
                rank_logl(iproc), rank_aic(iproc), rank_bic(iproc)
        else
            write(*, '(A8,ES12.3,F8.4,F9.4,F9.4,F7.2,F9.2,F12.1,F12.1,F12.1,I6,3I6)') &
                trim(proc_names(iproc)), omegas(iproc), alphas(iproc), gamma_ps(iproc), &
                betas(iproc), nus(iproc), vol_anns(iproc), &
                logls(iproc), aics(iproc), bics(iproc), niters(iproc), &
                rank_logl(iproc), rank_aic(iproc), rank_bic(iproc)
        end if
        if (.not. conv(iproc)) write(*, '(4X,A)') "WARNING: did not converge"
    end do

    write(*,*)

end do

write(*, '(A)') "par3 = theta (NAGARCH leverage shift), gamma (GJR/EGARCH asymmetry)"
write(*, '(A,F0.0,A)') "vol_ann% = sqrt(", trading_days, "*h_unc)*100 for GARCH/NAGARCH/GJR;"
write(*, '(A,F0.0,A)') "          sqrt(", trading_days, "*exp(omega/(1-beta)))*100 for EGARCH"

call cpu_time(t_end)
write(*, '(/,A,F8.3,A)') "elapsed time: ", t_end - t_start, " s"

end program xfit_garch_returns

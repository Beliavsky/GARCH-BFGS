program xfit_garch_dist_returns
! Read price data from a CSV and fit all combinations of
!   4 GARCH processes × selected noise distributions
! to the log-returns of each column.
! Reports parameter estimates, annualised unconditional vol,
! log-likelihood, AIC, BIC, and model ranks per asset.

use kind_mod,       only: dp
use csv_mod,        only: read_price_csv
use garch_flex_mod, only: flex_set_data, flex_set_types, flex_np, flex_obj, &
                          proc_garch, proc_nagarch, proc_gjr, proc_egarch, &
                           proc_names, &
                           dist_normal, dist_t, dist_sech, dist_ged, &
                           dist_laplace, dist_logistic, dist_nig
use garch_module,   only: garch_inv_transform,   garch_transform
use nagarch_module, only: nagarch_inv_transform, nagarch_transform
use gjr_module,     only: gjr_inv_transform,     gjr_transform
use egarch_module,  only: egarch_inv_transform,  egarch_transform
use bfgs_module,    only: bfgs_minimize
use stats_mod,      only: mean, sd
use rank_mod,       only: rank_desc, rank_asc
implicit none

integer , parameter :: n_proc       = proc_egarch ! 4 process types
integer , parameter :: n_dist_id    = dist_nig    ! total number of distribution IDs
real(dp), parameter :: trading_days = 252.0_dp
integer , parameter :: max_iter     = 100         ! maximum BFGS iterations per fit
real(dp), parameter :: gtol         = 1.0e-7_dp   ! BFGS convergence tolerance on gradient norm
integer , parameter :: nret         = 10**6       ! cap on observations used
logical , parameter :: log_ret      = .true.      ! .true. = log returns, .false. = simple returns
integer , parameter :: len_name     = 8
! distribution tables (row index = dist ID 1..n_dist_id)
character(len=len_name), parameter :: all_dist_names(n_dist_id) = [character(len=len_name) :: &
    "Normal", "t", "Sech", "GED", "Laplace", "Logistic", "NIG"]
! which distributions have a single shape parameter
logical, parameter :: has_shape(n_dist_id) = &
    [.false., .true., .false., .true., .false., .false., .true.]
! unconstrained starting value for shape parameter where applicable (ignored otherwise)
!   t:   nu ≈ 8    → p = -log(98/6 - 1)      ≈ -2.729
!   GED: nu = 1.5  → p = log(1.5)             ≈  0.405
!   NIG: alp = 2.0 → p = logit((2.0-0.1)/19.9) ≈ -2.249
real(dp), parameter :: dist_p0(n_dist_id) = &
    [0.0_dp, -2.729_dp, 0.0_dp, 0.405_dp, 0.0_dp, 0.0_dp, -2.249_dp]

! selected distributions to fit (edit this list to restrict the search)
integer, parameter :: sel_dists(*) = [dist_normal, dist_t, dist_sech, dist_ged, dist_laplace, dist_logistic, dist_nig]
integer, parameter :: n_dist = size(sel_dists), n_combo = n_proc * n_dist
! data
character(len=*), parameter :: prices_file = "vix_spy.csv"
integer,           allocatable :: dates(:)
character(len=32), allocatable :: col_names(:)
real(dp),          allocatable :: prices(:,:)
real(dp),          allocatable :: raw_ret(:)      ! raw returns before demeaning
real(dp),          allocatable :: ret(:)         ! demeaned returns
integer :: nobs, nprices, nall, ncols, icol, i1_price

real(dp) :: ret_mean, ret_std

! optimisation
real(dp), allocatable :: p(:), p0(:)
real(dp) :: fopt
integer  :: np, niter, iproc, id, idist, im
logical  :: converged
real(dp) :: t_start, t_end

! per-combo results (index im = (iproc-1)*n_dist + id)
real(dp), allocatable :: omegas(:), alphas(:), gamma_ps(:), betas(:)
real(dp), allocatable :: shapes(:), vol_anns(:), logls(:), aics(:), bics(:)
integer,  allocatable :: niters(:)
logical,  allocatable :: conv(:)

! rank arrays
integer, allocatable :: rank_logl(:), rank_aic(:), rank_bic(:)

! working scalars
real(dp) :: omega, alpha, gamma_p, beta, theta, shape_val, h_unc, vol_ann, logl, aic, bic


! single print format: gamma_p shown for all (0 for GARCH), shape shown for all (0 if N/A)
character(len=*), parameter :: row_fmt = &
    "(A8,A9,ES12.3,F8.4,F9.4,F9.4,F7.2,F9.2,F12.1,F12.1,F12.1,I6,3I6)"

call cpu_time(t_start)

! ── initialise ───────────────────────────────────────────────────────────────

allocate(omegas(n_combo), alphas(n_combo), gamma_ps(n_combo), betas(n_combo))
allocate(shapes(n_combo), vol_anns(n_combo), logls(n_combo), aics(n_combo), bics(n_combo))
allocate(niters(n_combo), conv(n_combo))
allocate(rank_logl(n_combo), rank_aic(n_combo), rank_bic(n_combo))


! ── read data ────────────────────────────────────────────────────────────────
call read_price_csv(prices_file, dates, col_names, prices)
nprices = size(prices, 1)
ncols   = size(col_names)
nall    = nprices - 1
nobs    = min(nret, nall)
i1_price = nprices-nobs
allocate(raw_ret(nobs), ret(nobs))

write(*, '(A,I0,A,I0,A,A,/)') "Using last ", nobs, " of ", nall, " observations, ", &
    merge("log returns   ", "simple returns", log_ret)

! ── loop over assets ─────────────────────────────────────────────────────────
do icol = 1, ncols

    if (log_ret) then
        raw_ret = log(prices(i1_price+1:nprices, icol) / prices(i1_price:nprices-1, icol))
    else
        raw_ret = prices(i1_price+1:nprices, icol) / prices(i1_price:nprices-1, icol) - 1.0_dp
    end if
    ret_mean = mean(raw_ret)
    ret      = raw_ret - ret_mean
    ret_std  = sd(ret)
    write(*, '(A,A)') "Asset: ", trim(col_names(icol))
    write(*, '(A,F8.4,A,F8.4,A)') &
        "  mean (bps): ", ret_mean*1.0e4_dp, "   emp vol_ann%: ", ret_std*sqrt(trading_days)*100.0_dp, "%"
    write(*, '(A8,A9,A12,A8,A9,A9,A7,A9,A12,A12,A12,A6,A8,A6,A6)') &
        "Process", "Dist", "omega", "alpha", "par3", "beta", "shape", &
        "vol_ann%", "logL", "AIC", "BIC", "iter", "rklogL", "rkAIC", "rkBIC"
    write(*, '(A)') repeat("-", 131)

    ! ── fit all process × distribution combinations ──────────────────────────
    do iproc = proc_garch, proc_egarch
        do id = 1, n_dist
            idist = sel_dists(id)
            im    = (iproc-1)*n_dist + id

            call flex_set_data(ret, nobs)
            call flex_set_types(iproc, idist)
            np = flex_np()
            allocate(p(np), p0(np))

            select case (iproc)
            case (proc_garch)
                call garch_inv_transform(1.0e-5_dp, 0.08_dp, 0.88_dp, p0(1:3))

            case (proc_nagarch)
                call nagarch_inv_transform(1.0e-5_dp, 0.06_dp, 0.88_dp, 0.5_dp, p0(1:4))

            case (proc_gjr)
                call gjr_inv_transform(1.0e-5_dp, 0.04_dp, 0.08_dp, 0.88_dp, p0(1:4))

            case (proc_egarch)
                call egarch_inv_transform(log(ret_std**2)*0.03_dp, 0.10_dp, -0.10_dp, 0.97_dp, p0(1:4))
            end select

            if (has_shape(idist)) p0(np) = dist_p0(idist)

            p = p0
            call bfgs_minimize(flex_obj, p, np, max_iter, gtol, fopt, niter, converged)

            ! extract shape parameter (default 0 for distributions without one)
            shape_val = 0.0_dp
            if (has_shape(idist)) then
                select case (idist)
                case (dist_t)
                    shape_val = 2.0_dp + 98.0_dp / (1.0_dp + exp(-p(np)))
                case (dist_ged)
                    shape_val = exp(p(np))
                case (dist_nig)
                    shape_val = 0.1_dp + 19.9_dp / (1.0_dp + exp(-p(np)))
                end select
            end if

            ! extract process parameters
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

            omegas(im)   = omega
            alphas(im)   = alpha
            gamma_ps(im) = gamma_p
            betas(im)    = beta
            shapes(im)   = shape_val
            vol_anns(im) = vol_ann
            logls(im)    = logl
            aics(im)     = aic
            bics(im)     = bic
            niters(im)   = niter
            conv(im)     = converged

            deallocate(p, p0)

        end do
    end do

    call rank_desc(logls, rank_logl)
    call rank_asc(aics,   rank_aic)
    call rank_asc(bics,   rank_bic)

    ! ── print results ────────────────────────────────────────────────────────
    do iproc = proc_garch, proc_egarch
        do id = 1, n_dist
            idist = sel_dists(id)
            im    = (iproc-1)*n_dist + id
            write(*, row_fmt) &
                proc_names(iproc), all_dist_names(idist), &
                omegas(im), alphas(im), gamma_ps(im), betas(im), shapes(im), &
                vol_anns(im), logls(im), aics(im), bics(im), niters(im), &
                rank_logl(im), rank_aic(im), rank_bic(im)
            if (.not. conv(im)) write(*, '(4X,A)') "WARNING: did not converge"
        end do
    end do

    write(*, *)

end do

write(*, '(A)') "par3   = theta (NAGARCH leverage shift), gamma (GJR/EGARCH asymmetry); 0 for GARCH"
write(*, '(A)') "shape  = nu (t,GED), alpha (NIG); 0 = not applicable"
write(*, '(A,F0.0,A)') "vol_ann% = sqrt(", trading_days, "*h_unc)*100 for GARCH/NAGARCH/GJR;"
write(*, '(A,F0.0,A)') "          sqrt(", trading_days, "*exp(omega/(1-beta)))*100 for EGARCH"

call cpu_time(t_end)
write(*, '(/,A,F0.3,A)') "elapsed time: ", t_end - t_start, " s"

end program xfit_garch_dist_returns

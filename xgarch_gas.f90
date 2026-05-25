program xgarch_gas
! Fits GARCH and GAS volatility models to log-returns of each asset and
! reports joint rankings across all models per asset.
!
! GARCH: 4 processes (GARCH, NAGARCH, GJR, EGARCH) × sel_dists
! GAS:   2 processes (GAS, AGAS)                   × gas_sel_dists
!
! To change which distributions GARCH fits, edit sel_dists.
! To change which distributions GAS fits, edit gas_sel_dists
!   (currently restricted to Normal and t, which gas_mod supports).

use kind_mod,       only: dp
use csv_mod,        only: read_price_csv, print_price_sample_info
! dist_sech, dist_ged, dist_laplace, dist_logistic, dist_nig are also available
! from garch_flex_mod — add them here and to sel_dists to fit more distributions
use garch_flex_mod, only: flex_set_data, flex_set_types, flex_np, flex_obj, &
                           proc_garch, proc_nagarch, proc_gjr, proc_egarch, &
                           garch_proc_names => proc_names, &
                           dist_normal, dist_t, dist_ged, dist_nig
use garch_mod,   only: garch_inv_transform,   garch_transform
use nagarch_mod, only: nagarch_inv_transform, nagarch_transform
use gjr_mod,     only: gjr_inv_transform,     gjr_transform
use egarch_mod,  only: egarch_inv_transform,  egarch_transform
use gas_mod,        only: gas_set_data, gas_set_types, gas_np, gas_obj, &
                           gas_sym_inv_transform, gas_asym_inv_transform, gas_transform, &
                           proc_gas, proc_agas, &
                           gas_proc_names => proc_names, &
                           n_gas_proc => n_proc
use bfgs_mod,    only: bfgs_minimize
use stats_mod,      only: mean, sd
use rank_mod,       only: rank_desc, rank_asc
implicit none

integer,  parameter :: max_iter     = 100          ! maximum BFGS iterations per fit
real(dp), parameter :: gtol         = 1.0e-7_dp    ! BFGS convergence tolerance on gradient norm
integer,  parameter :: nret         = 10**6         ! cap on observations used
logical,  parameter :: log_ret      = .true.        ! .true. = log returns, .false. = simple returns
real(dp), parameter :: trading_days = 252.0_dp
integer,  parameter :: len_name     = 8
integer,  parameter :: n_dist_id    = dist_nig      ! total number of distribution IDs (7)

! distribution tables (indexed by dist ID 1..n_dist_id)
character(len=len_name), parameter :: all_dist_names(n_dist_id) = [character(len=len_name) :: &
    "Normal", "t", "Sech", "GED", "Laplace", "Logistic", "NIG"]
logical,  parameter :: has_shape(n_dist_id) = &
    [.false., .true., .false., .true., .false., .false., .true.]
! unconstrained starting value for shape parameter (ignored if has_shape is .false.)
real(dp), parameter :: dist_p0(n_dist_id) = &
    [0.0_dp, -2.729_dp, 0.0_dp, 0.405_dp, 0.0_dp, 0.0_dp, -2.249_dp]

! distributions to fit in GARCH models (edit to change)
integer, parameter :: sel_dists(*)  = [dist_normal, dist_t]
integer, parameter :: n_dist        = size(sel_dists)

! distributions to fit in GAS models (edit to change; must be supported by gas_mod)
integer, parameter :: gas_sel_dists(*) = [dist_normal, dist_t]
integer, parameter :: n_gas_dist       = size(gas_sel_dists)

integer, parameter :: n_garch_proc  = proc_egarch               ! 4
integer, parameter :: n_garch_combo = n_garch_proc * n_dist
integer, parameter :: n_gas_combo   = n_gas_proc   * n_gas_dist
integer, parameter :: n_combo       = n_garch_combo + n_gas_combo

character(len=*), parameter :: prices_file = "vix_spy.csv"
character(len=*), parameter :: row_fmt = &
    "(A8,A9,ES12.3,F8.4,F9.4,F9.4,F7.2,F9.2,F12.1,F12.1,F12.1,I6,3I6)"

! data
integer,           allocatable :: dates(:)
character(len=32), allocatable :: col_names(:)
real(dp),          allocatable :: prices(:,:)
real(dp),          allocatable :: raw_ret(:)  ! raw returns before demeaning
real(dp),          allocatable :: ret(:)      ! demeaned returns
integer :: nobs, nprices, nall, ncols, icol, i1_price

real(dp) :: ret_mean, ret_std

! optimisation
real(dp), allocatable :: p(:), p0(:)
real(dp) :: fopt
integer  :: np, niter, iproc, id, idist, im
logical  :: converged
real(dp) :: t_start, t_end

! per-combo results
real(dp), allocatable :: omegas(:), alphas(:), par3s(:), betas(:)
real(dp), allocatable :: shapes(:), vol_anns(:), logls(:), aics(:), bics(:)
integer,  allocatable :: niters(:)
logical,  allocatable :: conv(:)
integer,  allocatable :: rank_logl(:), rank_aic(:), rank_bic(:)

! working scalars
real(dp) :: omega, alpha, par3, beta, theta, shape_val, h_unc, vol_ann, logl, aic, bic

call cpu_time(t_start)

allocate(omegas(n_combo), alphas(n_combo), par3s(n_combo), betas(n_combo))
allocate(shapes(n_combo), vol_anns(n_combo), logls(n_combo), aics(n_combo), bics(n_combo))
allocate(niters(n_combo), conv(n_combo))
allocate(rank_logl(n_combo), rank_aic(n_combo), rank_bic(n_combo))

call read_price_csv(prices_file, dates, col_names, prices)
nprices  = size(prices, 1)
ncols    = size(col_names)
nall     = nprices - 1
nobs     = min(nret, nall)
i1_price = nprices - nobs
allocate(raw_ret(nobs), ret(nobs))

call print_price_sample_info(prices_file, dates, ncols, nobs, &
    merge("demeaned log returns   ", "demeaned simple returns", log_ret))
write(*,*)

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

    ! ── fit GARCH models ──────────────────────────────────────────────────────
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

            h_unc = 0.0_dp
            select case (iproc)
            case (proc_garch)
                call garch_transform(p(1:3), omega, alpha, beta)
                par3  = 0.0_dp
                h_unc = omega / (1.0_dp - alpha - beta)
            case (proc_nagarch)
                call nagarch_transform(p(1:4), omega, alpha, beta, theta)
                par3  = theta
                h_unc = omega / (1.0_dp - alpha*(1.0_dp + theta**2) - beta)
            case (proc_gjr)
                call gjr_transform(p(1:4), omega, alpha, par3, beta)
                h_unc = omega / (1.0_dp - alpha - 0.5_dp*par3 - beta)
            case (proc_egarch)
                call egarch_transform(p(1:4), omega, alpha, par3, beta)
                h_unc = exp(omega / (1.0_dp - beta))
            end select

            vol_ann      = sqrt(trading_days * h_unc) * 100.0_dp
            logl         = -nobs * fopt
            aic          = 2.0_dp*np             - 2.0_dp*logl
            bic          = np*log(real(nobs,dp)) - 2.0_dp*logl
            omegas(im)   = omega
            alphas(im)   = alpha
            par3s(im)    = par3
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

    ! ── fit GAS models ────────────────────────────────────────────────────────
    do iproc = proc_gas, proc_agas
        do id = 1, n_gas_dist
            idist = gas_sel_dists(id)
            im    = n_garch_combo + (iproc-1)*n_gas_dist + id

            call gas_set_data(ret, nobs)
            call gas_set_types(iproc, idist)
            np = gas_np()
            allocate(p(np), p0(np))

            select case (iproc)
            case (proc_gas)
                call gas_sym_inv_transform(log(ret_std**2)*0.03_dp, 0.1_dp, 0.97_dp, p0(1:3))
            case (proc_agas)
                call gas_asym_inv_transform(log(ret_std**2)*0.03_dp, 0.1_dp, 0.0_dp, 0.97_dp, p0(1:4))
            end select
            if (has_shape(idist)) p0(np) = dist_p0(idist)

            p = p0
            call bfgs_minimize(gas_obj, p, np, max_iter, gtol, fopt, niter, converged)

            call gas_transform(p, omega, alpha, par3, beta)
            shape_val = 0.0_dp
            if (has_shape(idist)) shape_val = 2.0_dp + 98.0_dp / (1.0_dp + exp(-p(np)))

            h_unc        = exp(omega / (1.0_dp - beta))
            vol_ann      = sqrt(trading_days * h_unc) * 100.0_dp
            logl         = -nobs * fopt
            aic          = 2.0_dp*np             - 2.0_dp*logl
            bic          = np*log(real(nobs,dp)) - 2.0_dp*logl
            omegas(im)   = omega
            alphas(im)   = alpha
            par3s(im)    = par3
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

    ! ── joint rankings ────────────────────────────────────────────────────────
    call rank_desc(logls, rank_logl)
    call rank_asc(aics,   rank_aic)
    call rank_asc(bics,   rank_bic)

    ! ── print GARCH block ─────────────────────────────────────────────────────
    do iproc = proc_garch, proc_egarch
        do id = 1, n_dist
            idist = sel_dists(id)
            im    = (iproc-1)*n_dist + id
            write(*, row_fmt) &
                garch_proc_names(iproc), all_dist_names(idist), &
                omegas(im), alphas(im), par3s(im), betas(im), shapes(im), &
                vol_anns(im), logls(im), aics(im), bics(im), niters(im), &
                rank_logl(im), rank_aic(im), rank_bic(im)
            if (.not. conv(im)) write(*, '(4X,A)') "WARNING: did not converge"
        end do
    end do

    ! ── separator then GAS block ──────────────────────────────────────────────
    write(*, '(A)') repeat("-", 131)

    do iproc = proc_gas, proc_agas
        do id = 1, n_gas_dist
            idist = gas_sel_dists(id)
            im    = n_garch_combo + (iproc-1)*n_gas_dist + id
            write(*, row_fmt) &
                gas_proc_names(iproc), all_dist_names(idist), &
                omegas(im), alphas(im), par3s(im), betas(im), shapes(im), &
                vol_anns(im), logls(im), aics(im), bics(im), niters(im), &
                rank_logl(im), rank_aic(im), rank_bic(im)
            if (.not. conv(im)) write(*, '(4X,A)') "WARNING: did not converge"
        end do
    end do

    write(*, *)

end do

write(*, '(A)') "par3     = theta (NAGARCH leverage), gamma (GJR/EGARCH/AGAS asymmetry); 0 for GARCH/GAS"
write(*, '(A)') "shape    = nu (t,GED), alpha (NIG); 0 = not applicable"
write(*, '(A,F0.0,A)') "vol_ann% = sqrt(", trading_days, "*h_unc)*100"
write(*, '(A)') "         GARCH/NAGARCH/GJR: h_unc = omega / (1 - alpha*(1+par3^2) - beta)"
write(*, '(A)') "         EGARCH/GAS/AGAS:   h_unc = exp(omega / (1 - beta))"

call cpu_time(t_end)
write(*, '(/,A,F0.3,A)') "elapsed time: ", t_end - t_start, " s"

end program xgarch_gas

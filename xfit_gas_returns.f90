program xfit_gas_returns
! Read price data from a CSV and fit all combinations of
!   2 GAS processes × 2 noise distributions
! to the log-returns of each column.
! Reports parameter estimates, annualised unconditional vol,
! log-likelihood, AIC, BIC, and model ranks per asset.

use kind_mod,    only: dp
use csv_mod,     only: read_price_csv
use gas_mod,     only: gas_set_data, gas_set_types, gas_np, gas_obj, &
                        gas_sym_inv_transform, gas_asym_inv_transform, gas_transform, &
                        proc_gas, proc_agas, dist_t, &
                        n_proc, n_dist, proc_names, dist_names, has_shape
use bfgs_module, only: bfgs_minimize
use stats_mod,   only: mean, sd
use rank_mod,    only: rank_desc, rank_asc
implicit none

integer,  parameter :: max_iter     = 100          ! maximum BFGS iterations per fit
real(dp), parameter :: gtol         = 1.0e-7_dp    ! BFGS convergence tolerance on gradient norm
integer,  parameter :: nret         = 10**6         ! cap on observations used
logical,  parameter :: log_ret      = .true.        ! .true. = log returns, .false. = simple returns
real(dp), parameter :: trading_days = 252.0_dp
integer,  parameter :: n_combo      = n_proc * n_dist
real(dp), parameter :: nu_p0        = -2.729_dp    ! t-dist start: nu ~ 8

character(len=*), parameter :: prices_file = "vix_spy.csv"

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

! per-combo results (index im = (iproc-1)*n_dist + id)
real(dp), allocatable :: omegas(:), alphas(:), gammas(:), betas(:)
real(dp), allocatable :: shapes(:), vol_anns(:), logls(:), aics(:), bics(:)
integer,  allocatable :: niters(:)
logical,  allocatable :: conv(:)

! rank arrays
integer, allocatable :: rank_logl(:), rank_aic(:), rank_bic(:)

! working scalars
real(dp) :: omega, alpha, gamma, beta, h_unc, vol_ann, logl, aic, bic, shape_val

character(len=*), parameter :: row_fmt = &
    "(A8,A9,ES12.3,F8.4,F9.4,F9.4,F7.2,F9.2,F12.1,F12.1,F12.1,I6,3I6)"

call cpu_time(t_start)

allocate(omegas(n_combo), alphas(n_combo), gammas(n_combo), betas(n_combo))
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

write(*, '(A,I0,A,I0,A,A,/)') "Using last ", nobs, " of ", nall, " observations, ", &
    merge("log returns   ", "simple returns", log_ret)

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
        "Process", "Dist", "omega", "alpha", "gamma", "beta", "shape", &
        "vol_ann%", "logL", "AIC", "BIC", "iter", "rklogL", "rkAIC", "rkBIC"
    write(*, '(A)') repeat("-", 131)

    do iproc = proc_gas, proc_agas
        do id = 1, n_dist
            idist = id
            im    = (iproc-1)*n_dist + id

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
            if (has_shape(idist)) p0(np) = nu_p0

            p = p0
            call bfgs_minimize(gas_obj, p, np, max_iter, gtol, fopt, niter, converged)

            call gas_transform(p, omega, alpha, gamma, beta)

            shape_val = 0.0_dp
            if (idist == dist_t) shape_val = 2.0_dp + 98.0_dp / (1.0_dp + exp(-p(np)))

            h_unc   = exp(omega / (1.0_dp - beta))
            vol_ann = sqrt(trading_days * h_unc) * 100.0_dp
            logl    = -nobs * fopt
            aic     = 2.0_dp*np             - 2.0_dp*logl
            bic     = np*log(real(nobs,dp)) - 2.0_dp*logl

            omegas(im)   = omega
            alphas(im)   = alpha
            gammas(im)   = gamma
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

    do iproc = proc_gas, proc_agas
        do id = 1, n_dist
            im = (iproc-1)*n_dist + id
            write(*, row_fmt) &
                proc_names(iproc), dist_names(id), &
                omegas(im), alphas(im), gammas(im), betas(im), shapes(im), &
                vol_anns(im), logls(im), aics(im), bics(im), niters(im), &
                rank_logl(im), rank_aic(im), rank_bic(im)
            if (.not. conv(im)) write(*, '(4X,A)') "WARNING: did not converge"
        end do
    end do

    write(*, *)

end do

write(*, '(A)') "gamma    = asymmetry coefficient (GJR-style score weight for negative shocks; 0 for GAS)"
write(*, '(A)') "shape    = nu (t-dist degrees of freedom); 0 = not applicable"
write(*, '(A,F0.0,A)') "vol_ann% = sqrt(", trading_days, "*exp(omega/(1-beta)))*100"

call cpu_time(t_end)
write(*, '(/,A,F0.3,A)') "elapsed time: ", t_end - t_start, " s"

end program xfit_gas_returns

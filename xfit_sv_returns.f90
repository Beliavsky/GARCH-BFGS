program xfit_sv_returns
! Read price data from a CSV and fit all four SV variants to log-returns:
!   SV-N, SV-lev-N, SV-t, SV-lev-t
! via QML (Kalman filter).  Reports parameter estimates, annualised
! unconditional vol, log-likelihood, AIC, BIC, and model ranks.
! Note: QML log-likelihoods are not comparable to GARCH/GAS values.

use kind_mod,    only: dp
use csv_mod,     only: read_price_csv
use sv_mod,      only: sv_set_data, sv_set_types, sv_np, sv_obj, sv_skew_kurt, &
                        sv_sym_inv_transform, sv_lev_inv_transform, &
                        sv_t_inv_transform,   sv_lev_t_inv_transform, &
                        sv_transform, &
                        proc_sv, proc_sv_lev, dist_normal, dist_t, &
                        n_proc, n_dist, model_names
use bfgs_module, only: bfgs_minimize
use stats_mod,   only: mean, sd
use rank_mod,    only: rank_desc, rank_asc
implicit none

integer, parameter :: nmod = n_proc * n_dist   ! 4 models

! flat model index: imod = (iproc-1)*n_dist + idist
! imod: 1=(SV-N), 2=(SV-t), 3=(SV-lev-N), 4=(SV-lev-t)

! constants
real(dp), parameter :: trading_days = 252.0_dp
real(dp), parameter :: nu0          = 8.0_dp    ! t starting value
integer,  parameter :: max_iter     = 500
real(dp), parameter :: gtol         = 1.0e-4_dp
integer,  parameter :: nret         = 10**4

! data
integer,           allocatable :: dates(:)
character(len=32), allocatable :: col_names(:)
real(dp),          allocatable :: prices(:,:)
real(dp),          allocatable :: raw_ret(:), ret(:)
integer :: nobs, nprices, nall, ncols, icol

real(dp) :: ret_mean, ret_std

! optimisation
real(dp), allocatable :: p(:), p0(:)
real(dp) :: fopt
integer  :: np, niter, iproc, idist, imod
logical  :: converged
real(dp) :: t_start, t_end

! per-model results
real(dp) :: mus(nmod), phis(nmod), ses(nmod), rhos(nmod), nus(nmod)
real(dp) :: vol_anns(nmod), logls(nmod), aics(nmod), bics(nmod)
real(dp) :: skews(nmod), kurts(nmod)
integer  :: niters(nmod)
logical  :: conv(nmod)

! working scalars
real(dp) :: mu, phi, sigma_eta, rho, nu, h_unc, vol_ann, logl, aic, bic
character(len=*), parameter :: prices_file = "vix_spy.csv"
integer :: rank_logl(nmod), rank_aic(nmod), rank_bic(nmod)

! format strings: model(10), mu(9), phi(8), sig_eta(9), rho(8|8X), nu(7|7X), vol%(9), logL(12), AIC(12), BIC(12), iter(6), skew(8), ekurt(8), ranks(3x6)
character(len=*), parameter :: fmt_nn  = "(A10,F9.4,F8.4,F9.4,8X,7X,F9.2,F12.1,F12.1,F12.1,I6,F8.3,F8.3,3I6)"
character(len=*), parameter :: fmt_nt  = "(A10,F9.4,F8.4,F9.4,8X,F7.2,F9.2,F12.1,F12.1,F12.1,I6,F8.3,F8.3,3I6)"
character(len=*), parameter :: fmt_ln  = "(A10,F9.4,F8.4,F9.4,F8.4,7X,F9.2,F12.1,F12.1,F12.1,I6,F8.3,F8.3,3I6)"
character(len=*), parameter :: fmt_lt  = "(A10,F9.4,F8.4,F9.4,F8.4,F7.2,F9.2,F12.1,F12.1,F12.1,I6,F8.3,F8.3,3I6)"

call read_price_csv(prices_file, dates, col_names, prices)
nprices = size(prices, 1)
ncols   = size(col_names)
nall    = nprices - 1
nobs    = min(nret, nall)
allocate(raw_ret(nobs), ret(nobs))

write(*, '(A,I0,A,I0,A)') "Using last ", nobs, " of ", nall, " observations"
write(*,*)

call cpu_time(t_start)

do icol = 1, ncols

    raw_ret  = log(prices(nall-nobs+2:nprices, icol) / prices(nall-nobs+1:nprices-1, icol))
    ret_mean = mean(raw_ret)
    ret      = raw_ret - ret_mean
    ret_std  = sd(ret)

    write(*, '(A,A)') "Asset: ", trim(col_names(icol))
    write(*, '(A,F8.4,A,F8.4,A)') &
        "  mean (bps): ", ret_mean*1.0e4_dp, "   emp vol_ann%: ", ret_std*sqrt(trading_days)*100.0_dp, "%"
    write(*, '(A10,A9,A8,A9,A8,A7,A9,A12,A12,A12,A6,A8,A8,A8,A6,A6)') &
        "Model", "mu", "phi", "sig_eta", "rho", "nu", &
        "vol_ann%", "logL", "AIC", "BIC", "iter", "skew", "ekurt", "rklogL", "rkAIC", "rkBIC"
    write(*, '(A)') repeat("-", 136)

    call sv_set_data(ret, nobs)

    do iproc = proc_sv, proc_sv_lev
        do idist = dist_normal, dist_t

            imod = (iproc-1)*n_dist + idist

            call sv_set_types(iproc, idist)
            np = sv_np()
            allocate(p(np), p0(np))

            select case (iproc*10 + idist)
            case (proc_sv*10     + dist_normal)
                call sv_sym_inv_transform(log(ret_std**2), 0.97_dp, 0.20_dp, p0)
            case (proc_sv*10     + dist_t)
                call sv_t_inv_transform(log(ret_std**2), 0.97_dp, 0.20_dp, nu0, p0)
            case (proc_sv_lev*10 + dist_normal)
                call sv_lev_inv_transform(log(ret_std**2), 0.97_dp, 0.20_dp, -0.3_dp, p0)
            case (proc_sv_lev*10 + dist_t)
                call sv_lev_t_inv_transform(log(ret_std**2), 0.97_dp, 0.20_dp, -0.3_dp, nu0, p0)
            end select

            p = p0
            call bfgs_minimize(sv_obj, p, np, max_iter, gtol, fopt, niter, converged)
            call sv_transform(p, mu, phi, sigma_eta, rho, nu)

            h_unc   = exp(mu + sigma_eta**2 / (2.0_dp*(1.0_dp - phi**2)))
            vol_ann = sqrt(trading_days * h_unc) * 100.0_dp
            logl    = -nobs * fopt
            aic     = 2.0_dp*np             - 2.0_dp*logl
            bic     = np*log(real(nobs,dp)) - 2.0_dp*logl

            mus(imod)      = mu
            phis(imod)     = phi
            ses(imod)      = sigma_eta
            rhos(imod)     = rho
            nus(imod)      = nu
            vol_anns(imod) = vol_ann
            logls(imod)    = logl
            aics(imod)     = aic
            bics(imod)     = bic
            niters(imod)   = niter
            conv(imod)     = converged
            call sv_skew_kurt(p, np, skews(imod), kurts(imod))

            deallocate(p, p0)

        end do
    end do

    call rank_desc(logls, rank_logl)
    call rank_asc(aics,   rank_aic)
    call rank_asc(bics,   rank_bic)

    do iproc = proc_sv, proc_sv_lev
        do idist = dist_normal, dist_t
            imod = (iproc-1)*n_dist + idist
            select case (iproc*10 + idist)
            case (proc_sv*10     + dist_normal)
                write(*, fmt_nn) trim(model_names(iproc,idist)), &
                    mus(imod), phis(imod), ses(imod), &
                    vol_anns(imod), logls(imod), aics(imod), bics(imod), niters(imod), &
                    skews(imod), kurts(imod), rank_logl(imod), rank_aic(imod), rank_bic(imod)
            case (proc_sv*10     + dist_t)
                write(*, fmt_nt) trim(model_names(iproc,idist)), &
                    mus(imod), phis(imod), ses(imod), nus(imod), &
                    vol_anns(imod), logls(imod), aics(imod), bics(imod), niters(imod), &
                    skews(imod), kurts(imod), rank_logl(imod), rank_aic(imod), rank_bic(imod)
            case (proc_sv_lev*10 + dist_normal)
                write(*, fmt_ln) trim(model_names(iproc,idist)), &
                    mus(imod), phis(imod), ses(imod), rhos(imod), &
                    vol_anns(imod), logls(imod), aics(imod), bics(imod), niters(imod), &
                    skews(imod), kurts(imod), rank_logl(imod), rank_aic(imod), rank_bic(imod)
            case (proc_sv_lev*10 + dist_t)
                write(*, fmt_lt) trim(model_names(iproc,idist)), &
                    mus(imod), phis(imod), ses(imod), rhos(imod), nus(imod), &
                    vol_anns(imod), logls(imod), aics(imod), bics(imod), niters(imod), &
                    skews(imod), kurts(imod), rank_logl(imod), rank_aic(imod), rank_bic(imod)
            end select
            if (.not. conv(imod)) write(*, '(4X,A)') "WARNING: did not converge"
        end do
    end do

    write(*,*)

end do

write(*, '(A)') "vol_ann% = sqrt(252 * exp(mu + sig_eta^2 / (2*(1-phi^2)))) * 100"
write(*, '(A)') "Note: QML logL is not comparable to GARCH/GAS logL (different observation scale)"

call cpu_time(t_end)
write(*, '(/,A,F8.3,A)') "elapsed time: ", t_end - t_start, " s"

end program xfit_sv_returns

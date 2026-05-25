program xfit_sv_garch_returns
! Read price data from a CSV and fit, for each asset:
!   4 SV variants:   SV-N, SV-t, SV-lev-N, SV-lev-t  (QML)
!   8 GARCH variants: GARCH/NAGARCH/GJR/EGARCH x Normal/t
! Models are ranked within each class separately.
! WARNING: SV QML logL and GARCH logL are on different observation scales
!          and cannot be compared directly.

use kind_mod,         only: dp
use csv_mod,          only: read_price_csv, print_price_sample_info
use sv_mod,           only: sv_set_data, sv_set_types, sv_np, sv_obj, sv_pred_logl, sv_skew_kurt, &
                              sv_sym_inv_transform, sv_lev_inv_transform, &
                              sv_t_inv_transform,   sv_lev_t_inv_transform, &
                              sv_transform, &
                              proc_sv, proc_sv_lev, dist_normal, dist_t, &
                              sv_n_proc => n_proc, sv_n_dist => n_dist, &
                              sv_model_names => model_names
use garch_flex_mod,   only: flex_set_data, flex_set_types, flex_np, flex_obj, flex_skew_kurt, &
                              proc_nagarch, proc_anagarch, proc_cnagarch, dist_skew_t, &
                              anagarch_transform, anagarch_inv_transform, &
                              cnagarch_transform, cnagarch_inv_transform
use fgarch_mod,    only: fgarch_set_data, fgarch_set_dist, fgarch_np, fgarch_obj, &
                              fgarch_transform, fgarch_inv_transform, &
                              fgarch_vol_ann, fgarch_skew_kurt, &
                              fg_dist_normal, fg_dist_t
use nagarch_mod,   only: nagarch_inv_transform, nagarch_transform
use bfgs_mod,      only: bfgs_minimize
use stats_mod,        only: mean, sd
use rank_mod,         only: rank_desc, rank_asc
implicit none

integer, parameter :: sv_nmod         = sv_n_proc * sv_n_dist   ! 4
integer, parameter :: n_garch_dist    = 3
integer, parameter :: n_garch_procs   = 2   ! NAGARCH (squared) and ANAGARCH (abs-value)
integer, parameter :: garch_nmod      = n_garch_procs * n_garch_dist  ! 6
integer, parameter :: cn_nmod         = n_garch_dist                  ! Component NAGARCH × 3 dists
logical, parameter :: fit_fgarch      = .false.                         ! set .true. to enable FGARCH
integer, parameter :: fg_nmod         = 2                              ! FGARCH-N, FGARCH-t
integer, parameter :: nmod_total      = sv_nmod + garch_nmod + cn_nmod ! 4+6+3=13
integer, parameter :: garch_dists(n_garch_dist)       = [dist_normal, dist_t, dist_skew_t]
integer, parameter :: garch_procs(n_garch_procs)      = [proc_nagarch, proc_anagarch]

! Family GARCH model names
character(len=10), parameter :: fg_mnames(fg_nmod) = ["FGARCH-N  ", "FGARCH-t  "]
integer,           parameter :: fg_dists(fg_nmod)  = [fg_dist_normal, fg_dist_t]

! GARCH model names (n_garch_procs x {N, t, skew-t})
character(len=10), parameter :: garch_mnames(n_garch_procs, n_garch_dist) = reshape( &
    ["NAGARCH-N ", "ANAGARCH-N", &
     "NAGARCH-t ", "ANAGARCH-t", &
     "NAGARCH-sk", "ANAGARCH-s"], &
    [n_garch_procs, n_garch_dist])

! Component NAGARCH model names
character(len=10), parameter :: cn_mnames(cn_nmod) = &
    ["CNAGARCH-N", "CNAGARCH-t", "CNAGARCH-s"]

! constants
real(dp), parameter :: trading_days = 252.0_dp
real(dp), parameter :: nu0          = 8.0_dp     ! t starting value (nu ~ 8)
real(dp), parameter :: p_t0         = -2.729_dp  ! atanh-like transform for nu0
integer,  parameter :: max_iter     = 500
real(dp), parameter :: gtol         = 1.0e-4_dp
integer,  parameter :: nret         = 10**6

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
integer  :: i_gd   ! index into garch_dists(n_garch_dist)
integer  :: i_gp   ! index into garch_procs(n_garch_procs)
integer  :: i_cn   ! index into Component NAGARCH models

! SV per-model results
real(dp) :: sv_mus(sv_nmod), sv_phis(sv_nmod), sv_ses(sv_nmod)
real(dp) :: sv_rhos(sv_nmod), sv_nus(sv_nmod)
real(dp) :: sv_vol_anns(sv_nmod), sv_logls(sv_nmod), sv_aics(sv_nmod), sv_bics(sv_nmod)
real(dp) :: sv_pred_logls(sv_nmod)
real(dp) :: sv_skews(sv_nmod), sv_kurts(sv_nmod)
integer  :: sv_nps(sv_nmod)
integer  :: sv_niters(sv_nmod)
logical  :: sv_conv(sv_nmod)
integer  :: sv_rank_logl(sv_nmod), sv_rank_aic(sv_nmod), sv_rank_bic(sv_nmod)

! GARCH per-model results
real(dp) :: g_omegas(garch_nmod), g_alphas(garch_nmod), g_par3s(garch_nmod)
real(dp) :: g_betas(garch_nmod), g_nus(garch_nmod)
real(dp) :: g_skews(garch_nmod), g_kurts(garch_nmod), g_gammas(garch_nmod)
real(dp) :: g_vol_anns(garch_nmod), g_logls(garch_nmod), g_aics(garch_nmod), g_bics(garch_nmod)
integer  :: g_niters(garch_nmod)
logical  :: g_conv(garch_nmod)
integer  :: g_rank_logl(garch_nmod), g_rank_aic(garch_nmod), g_rank_bic(garch_nmod)

! Family GARCH per-model results
real(dp) :: fg_omegas(fg_nmod), fg_alphas(fg_nmod), fg_betas(fg_nmod)
real(dp) :: fg_lams(fg_nmod), fg_nus(fg_nmod), fg_bs(fg_nmod), fg_cs(fg_nmod)
real(dp) :: fg_nu_ts(fg_nmod)
real(dp) :: fg_vol_anns(fg_nmod), fg_logls(fg_nmod), fg_aics(fg_nmod), fg_bics(fg_nmod)
real(dp) :: fg_skews(fg_nmod), fg_kurts(fg_nmod)
integer  :: fg_nps(fg_nmod), fg_niters(fg_nmod)
logical  :: fg_conv(fg_nmod)
integer  :: fg_rank_logl(fg_nmod), fg_rank_aic(fg_nmod), fg_rank_bic(fg_nmod)

! Component NAGARCH per-model results
real(dp) :: cn_omegas(cn_nmod), cn_alphas(cn_nmod), cn_betas(cn_nmod)
real(dp) :: cn_rhos(cn_nmod),   cn_phis(cn_nmod),   cn_thetas(cn_nmod)
real(dp) :: cn_nus(cn_nmod),    cn_gammas(cn_nmod)
real(dp) :: cn_vol_anns(cn_nmod), cn_logls(cn_nmod), cn_aics(cn_nmod), cn_bics(cn_nmod)
real(dp) :: cn_skews(cn_nmod),  cn_kurts(cn_nmod)
integer  :: cn_niters(cn_nmod)
logical  :: cn_conv(cn_nmod)
integer  :: cn_rank_logl(cn_nmod), cn_rank_aic(cn_nmod), cn_rank_bic(cn_nmod)

! combined ranking
character(len=10) :: all_names(nmod_total)
real(dp) :: all_logls(nmod_total), all_aics(nmod_total), all_bics(nmod_total)
real(dp) :: all_skews(nmod_total), all_kurts(nmod_total)
integer  :: all_ranks(nmod_total), irank, j

! working scalars
real(dp) :: mu, phi, sigma_eta, rho, nu, gamma_val = 0.0_dp
real(dp) :: omega, alpha, g_par3 = 0.0_dp, beta, theta, h_unc, vol_ann, logl, aic, bic
real(dp) :: cn_rho, cn_phi, cn_D1, cn_Anum, cn_qbar, cn_hbar
real(dp) :: fg_omega, fg_alpha, fg_beta, fg_lam, fg_nu, fg_b, fg_c, fg_nu_t
real(dp) :: p0_fg(8)     ! full 8-element starting-value buffer for Family GARCH
integer  :: i_fg
character(len=*), parameter :: prices_file = "vix_spy.csv"

! ── SV format strings ────────────────────────────────────────────────────────
! columns: Model(10), mu(9), phi(8), sig_eta(9), rho(8|8X), nu(7|7X), vol%(9), logL(12), AIC(12), BIC(12), iter(6), skew(8), ekurt(8), 3x rank(6)
character(len=*), parameter :: sv_hdr  = "(A10,A9,A8,A9,A8,A7,A9,A12,A12,A12,A6,A8,A8,A6,A6,A6)"
character(len=*), parameter :: sv_nn   = "(A10,F9.4,F8.4,F9.4,8X,7X,F9.2,F12.1,F12.1,F12.1,I6,F8.3,F8.3,3I6)"
character(len=*), parameter :: sv_nt   = "(A10,F9.4,F8.4,F9.4,8X,F7.2,F9.2,F12.1,F12.1,F12.1,I6,F8.3,F8.3,3I6)"
character(len=*), parameter :: sv_ln   = "(A10,F9.4,F8.4,F9.4,F8.4,7X,F9.2,F12.1,F12.1,F12.1,I6,F8.3,F8.3,3I6)"
character(len=*), parameter :: sv_lt   = "(A10,F9.4,F8.4,F9.4,F8.4,F7.2,F9.2,F12.1,F12.1,F12.1,I6,F8.3,F8.3,3I6)"

! ── combined format strings ──────────────────────────────────────────────────
character(len=*), parameter :: comb_hdr = "(A6, 2X, A10, A12,  A12,  A12,  A8, A8)"
character(len=*), parameter :: comb_fmt = "(I6, 2X, A10, F12.1, F12.1, F12.1, F8.3, F8.3)"

! ── GARCH format strings ─────────────────────────────────────────────────────
! columns: Model(10), omega(12), alpha(8), par3(9|9X), beta(9), nu(7|7X), gamma(7|7X), vol%(9), logL(12), AIC(12), BIC(12), iter(6), skew(8), ekurt(8), 3x rank(6)
character(len=*), parameter :: g_hdr   = "(A10,A12,A8,A9,A9,A7,A7,A9,A12,A12,A12,A6,A8,A8,A6,A6,A6)"
character(len=*), parameter :: g_on    = "(A10,ES12.3,F8.4,F9.4,F9.4,7X,7X,F9.2,F12.1,F12.1,F12.1,I6,F8.3,F8.3,3I6)"
character(len=*), parameter :: g_ot    = "(A10,ES12.3,F8.4,F9.4,F9.4,F7.2,7X,F9.2,F12.1,F12.1,F12.1,I6,F8.3,F8.3,3I6)"
character(len=*), parameter :: g_osk   = "(A10,ES12.3,F8.4,F9.4,F9.4,F7.2,F7.3,F9.2,F12.1,F12.1,F12.1,I6,F8.3,F8.3,3I6)"

! ── Component NAGARCH format strings ─────────────────────────────────────────
! columns: Model(10), omega(12), alpha(8), beta(8), rho(8), phi(8), theta(8), nu(7|7X), gamma(7|7X), vol%(9), logL(12), AIC(12), BIC(12), iter(6), skew(8), ekurt(8), 3xrank(6)
character(len=*), parameter :: cn_hdr = &
    "(A10,A12,A8,A8,A8,A8,A8,A7,A7,A9,A12,A12,A12,A6,A8,A8,A6,A6,A6)"
character(len=*), parameter :: cn_n   = &
    "(A10,ES12.3,F8.4,F8.4,F8.4,F8.4,F8.4,7X,7X,F9.2,F12.1,F12.1,F12.1,I6,F8.3,F8.3,3I6)"
character(len=*), parameter :: cn_t   = &
    "(A10,ES12.3,F8.4,F8.4,F8.4,F8.4,F8.4,F7.2,7X,F9.2,F12.1,F12.1,F12.1,I6,F8.3,F8.3,3I6)"
character(len=*), parameter :: cn_sk  = &
    "(A10,ES12.3,F8.4,F8.4,F8.4,F8.4,F8.4,F7.2,F7.3,F9.2,F12.1,F12.1,F12.1,I6,F8.3,F8.3,3I6)"

! ── Family GARCH format strings ──────────────────────────────────────────────
! columns: Model(10), omega(10), alpha(7), beta(7), lam(7), nu(7), b(8), c(7), nu_t(7|7X), vol%(9), logL(12), AIC(12), BIC(12), iter(6), skew(8), ekurt(8), 3xrank(6)
character(len=*), parameter :: fg_hdr = &
    "(A10,A10,A7,A7,A7,A7,A8,A7,A7,A9,A12,A12,A12,A6,A8,A8,A6,A6,A6)"
character(len=*), parameter :: fg_fn  = &
    "(A10,ES10.3,F7.4,F7.4,F7.3,F7.3,F8.4,F7.4,7X,F9.2,F12.1,F12.1,F12.1,I6,F8.3,F8.3,3I6)"
character(len=*), parameter :: fg_ft  = &
    "(A10,ES10.3,F7.4,F7.4,F7.3,F7.3,F8.4,F7.4,F7.2,F9.2,F12.1,F12.1,F12.1,I6,F8.3,F8.3,3I6)"

call read_price_csv(prices_file, dates, col_names, prices)
nprices = size(prices, 1)
ncols   = size(col_names)
nall    = nprices - 1
nobs    = min(nret, nall)
allocate(raw_ret(nobs), ret(nobs))

call print_price_sample_info(prices_file, dates, ncols, nobs)
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

    ! ── SV models ────────────────────────────────────────────────────────────
    write(*, '(/,A)') "  SV models (QML):"
    write(*, sv_hdr) "Model", "mu", "phi", "sig_eta", "rho", "nu", &
        "vol_ann%", "logL", "AIC", "BIC", "iter", "skew", "ekurt", "rklogL", "rkAIC", "rkBIC"
    write(*, '(A)') "  " // repeat("-", 136)

    call sv_set_data(ret, nobs)

    do iproc = proc_sv, proc_sv_lev
        do idist = dist_normal, dist_t
            imod = (iproc-1)*sv_n_dist + idist
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

            sv_mus(imod)      = mu
            sv_phis(imod)     = phi
            sv_ses(imod)      = sigma_eta
            sv_rhos(imod)     = rho
            sv_nus(imod)      = nu
            sv_vol_anns(imod) = vol_ann
            sv_logls(imod)    = logl
            sv_aics(imod)     = aic
            sv_bics(imod)     = bic
            sv_niters(imod)   = niter
            sv_conv(imod)     = converged
            call sv_pred_logl(p, np, sv_pred_logls(imod))
            call sv_skew_kurt(p, np, sv_skews(imod), sv_kurts(imod))
            sv_nps(imod) = np

            deallocate(p, p0)
        end do
    end do

    call rank_desc(sv_logls, sv_rank_logl)
    call rank_asc(sv_aics,   sv_rank_aic)
    call rank_asc(sv_bics,   sv_rank_bic)

    do iproc = proc_sv, proc_sv_lev
        do idist = dist_normal, dist_t
            imod = (iproc-1)*sv_n_dist + idist
            select case (iproc*10 + idist)
            case (proc_sv*10     + dist_normal)
                write(*, sv_nn) trim(sv_model_names(iproc,idist)), &
                    sv_mus(imod), sv_phis(imod), sv_ses(imod), &
                    sv_vol_anns(imod), sv_logls(imod), sv_aics(imod), sv_bics(imod), &
                    sv_niters(imod), sv_skews(imod), sv_kurts(imod), &
                    sv_rank_logl(imod), sv_rank_aic(imod), sv_rank_bic(imod)
            case (proc_sv*10     + dist_t)
                write(*, sv_nt) trim(sv_model_names(iproc,idist)), &
                    sv_mus(imod), sv_phis(imod), sv_ses(imod), sv_nus(imod), &
                    sv_vol_anns(imod), sv_logls(imod), sv_aics(imod), sv_bics(imod), &
                    sv_niters(imod), sv_skews(imod), sv_kurts(imod), &
                    sv_rank_logl(imod), sv_rank_aic(imod), sv_rank_bic(imod)
            case (proc_sv_lev*10 + dist_normal)
                write(*, sv_ln) trim(sv_model_names(iproc,idist)), &
                    sv_mus(imod), sv_phis(imod), sv_ses(imod), sv_rhos(imod), &
                    sv_vol_anns(imod), sv_logls(imod), sv_aics(imod), sv_bics(imod), &
                    sv_niters(imod), sv_skews(imod), sv_kurts(imod), &
                    sv_rank_logl(imod), sv_rank_aic(imod), sv_rank_bic(imod)
            case (proc_sv_lev*10 + dist_t)
                write(*, sv_lt) trim(sv_model_names(iproc,idist)), &
                    sv_mus(imod), sv_phis(imod), sv_ses(imod), sv_rhos(imod), sv_nus(imod), &
                    sv_vol_anns(imod), sv_logls(imod), sv_aics(imod), sv_bics(imod), &
                    sv_niters(imod), sv_skews(imod), sv_kurts(imod), &
                    sv_rank_logl(imod), sv_rank_aic(imod), sv_rank_bic(imod)
            end select
            if (.not. sv_conv(imod)) write(*, '(4X,A)') "WARNING: did not converge"
        end do
    end do

    ! ── GARCH models ─────────────────────────────────────────────────────────
    write(*, '(/,A)') "  GARCH models:"
    write(*, g_hdr) "Model", "omega", "alpha", "par3", "beta", "nu", "gamma", &
        "vol_ann%", "logL", "AIC", "BIC", "iter", "skew", "ekurt", "rklogL", "rkAIC", "rkBIC"
    write(*, '(A)') "  " // repeat("-", 147)

    call flex_set_data(ret, nobs)

    do i_gp = 1, n_garch_procs
        iproc = garch_procs(i_gp)
        do i_gd = 1, n_garch_dist
            idist = garch_dists(i_gd)
            imod  = (i_gp - 1)*n_garch_dist + i_gd
            call flex_set_types(iproc, idist)
            np = flex_np()
            allocate(p(np), p0(np))

            ! starting values for structural parameters
            select case (iproc)
            case (proc_nagarch)
                call nagarch_inv_transform(1.0e-5_dp, 0.06_dp, 0.88_dp, 0.5_dp, p0(1:4))
            case (proc_anagarch)
                call anagarch_inv_transform(ret_std*0.07_dp, 0.10_dp, 0.85_dp, 0.0_dp, p0(1:4))
            end select
            if (idist == dist_t)      p0(np)   = p_t0
            if (idist == dist_skew_t) then
                p0(np-1) = p_t0
                p0(np)   = 0.0_dp
            end if

            p = p0
            call bfgs_minimize(flex_obj, p, np, max_iter, gtol, fopt, niter, converged)

            ! extract structural parameters
            select case (iproc)
            case (proc_nagarch)
                call nagarch_transform(p(1:4), omega, alpha, beta, theta)
                g_par3  = theta
                h_unc   = omega / (1.0_dp - alpha*(1.0_dp + theta**2) - beta)
                vol_ann = sqrt(trading_days * h_unc) * 100.0_dp
            case (proc_anagarch)
                ! abs-value: h_unc is unconditional sigma (not variance)
                call anagarch_transform(p(1:4), omega, alpha, beta, theta)
                g_par3  = theta
                h_unc   = omega / (1.0_dp - alpha - beta)
                vol_ann = h_unc * sqrt(trading_days) * 100.0_dp
            end select

            select case (idist)
            case (dist_normal)
                nu = 0.0_dp;  gamma_val = 0.0_dp
            case (dist_t)
                nu = 2.0_dp + 98.0_dp / (1.0_dp + exp(-p(np)));  gamma_val = 0.0_dp
            case (dist_skew_t)
                nu = 2.0_dp + 98.0_dp / (1.0_dp + exp(-p(np-1)));  gamma_val = exp(p(np))
            end select
            logl    = -nobs * fopt
            aic     = 2.0_dp*np             - 2.0_dp*logl
            bic     = np*log(real(nobs,dp)) - 2.0_dp*logl

            g_omegas(imod)   = omega
            g_alphas(imod)   = alpha
            g_par3s(imod)    = g_par3
            g_betas(imod)    = beta
            g_nus(imod)      = nu
            g_gammas(imod)   = gamma_val
            g_vol_anns(imod) = vol_ann
            g_logls(imod)    = logl
            g_aics(imod)     = aic
            g_bics(imod)     = bic
            g_niters(imod)   = niter
            g_conv(imod)     = converged
            call flex_skew_kurt(p, np, g_skews(imod), g_kurts(imod))

            deallocate(p, p0)
        end do
    end do

    call rank_desc(g_logls, g_rank_logl)
    call rank_asc(g_aics,   g_rank_aic)
    call rank_asc(g_bics,   g_rank_bic)

    ! Both NAGARCH and ANAGARCH have par3=theta; use g_on/g_ot/g_osk throughout
    do i_gp = 1, n_garch_procs
        iproc = garch_procs(i_gp)
        do i_gd = 1, n_garch_dist
            idist = garch_dists(i_gd)
            imod  = (i_gp - 1)*n_garch_dist + i_gd
            if (idist == dist_normal) then
                write(*, g_on) trim(garch_mnames(i_gp,i_gd)), &
                    g_omegas(imod), g_alphas(imod), g_par3s(imod), g_betas(imod), &
                    g_vol_anns(imod), g_logls(imod), g_aics(imod), g_bics(imod), &
                    g_niters(imod), g_skews(imod), g_kurts(imod), &
                    g_rank_logl(imod), g_rank_aic(imod), g_rank_bic(imod)
            else if (idist == dist_t) then
                write(*, g_ot) trim(garch_mnames(i_gp,i_gd)), &
                    g_omegas(imod), g_alphas(imod), g_par3s(imod), g_betas(imod), g_nus(imod), &
                    g_vol_anns(imod), g_logls(imod), g_aics(imod), g_bics(imod), &
                    g_niters(imod), g_skews(imod), g_kurts(imod), &
                    g_rank_logl(imod), g_rank_aic(imod), g_rank_bic(imod)
            else
                write(*, g_osk) trim(garch_mnames(i_gp,i_gd)), &
                    g_omegas(imod), g_alphas(imod), g_par3s(imod), g_betas(imod), &
                    g_nus(imod), g_gammas(imod), &
                    g_vol_anns(imod), g_logls(imod), g_aics(imod), g_bics(imod), &
                    g_niters(imod), g_skews(imod), g_kurts(imod), &
                    g_rank_logl(imod), g_rank_aic(imod), g_rank_bic(imod)
            end if
            if (.not. g_conv(imod)) write(*, '(4X,A)') "WARNING: did not converge"
        end do
    end do

    ! ── Component NAGARCH models ──────────────────────────────────────────────
    write(*, '(/,A)') "  Component NAGARCH models (Engle-Lee 1999):"
    write(*, cn_hdr) "Model", "omega", "alpha", "beta", "rho", "phi", "theta", &
        "nu", "gamma", "vol_ann%", "logL", "AIC", "BIC", "iter", "skew", "ekurt", &
        "rklogL", "rkAIC", "rkBIC"
    write(*, '(A)') "  " // repeat("-", 161)

    do i_cn = 1, cn_nmod
        idist = garch_dists(i_cn)
        call flex_set_types(proc_cnagarch, idist)
        np = flex_np()
        allocate(p(np), p0(np))

        call cnagarch_inv_transform(1.0e-6_dp, 0.06_dp, 0.88_dp, 0.99_dp, 0.01_dp, 0.5_dp, p0(1:6))
        if (idist == dist_t)      p0(np)   = p_t0
        if (idist == dist_skew_t) then
            p0(np-1) = p_t0
            p0(np)   = 0.0_dp
        end if

        p = p0
        call bfgs_minimize(flex_obj, p, np, max_iter, gtol, fopt, niter, converged)

        call cnagarch_transform(p(1:6), omega, alpha, beta, cn_rho, cn_phi, theta)

        ! Unconditional variance: solve coupled stationarity equations
        cn_D1   = max(1.0_dp - alpha*(1.0_dp + theta**2) - beta, 1.0e-8_dp)
        cn_Anum = (1.0_dp - alpha - beta) / cn_D1
        cn_qbar = omega / max((1.0_dp - cn_rho) - cn_phi * theta**2 * cn_Anum, 1.0e-8_dp)
        cn_hbar = cn_qbar * cn_Anum
        vol_ann = sqrt(trading_days * cn_hbar) * 100.0_dp

        select case (idist)
        case (dist_normal)
            nu = 0.0_dp;  gamma_val = 0.0_dp
        case (dist_t)
            nu = 2.0_dp + 98.0_dp / (1.0_dp + exp(-p(np)));  gamma_val = 0.0_dp
        case (dist_skew_t)
            nu = 2.0_dp + 98.0_dp / (1.0_dp + exp(-p(np-1)));  gamma_val = exp(p(np))
        end select
        logl = -nobs * fopt
        aic  = 2.0_dp*np             - 2.0_dp*logl
        bic  = np*log(real(nobs,dp)) - 2.0_dp*logl

        cn_omegas(i_cn)   = omega
        cn_alphas(i_cn)   = alpha
        cn_betas(i_cn)    = beta
        cn_rhos(i_cn)     = cn_rho
        cn_phis(i_cn)     = cn_phi
        cn_thetas(i_cn)   = theta
        cn_nus(i_cn)      = nu
        cn_gammas(i_cn)   = gamma_val
        cn_vol_anns(i_cn) = vol_ann
        cn_logls(i_cn)    = logl
        cn_aics(i_cn)     = aic
        cn_bics(i_cn)     = bic
        cn_niters(i_cn)   = niter
        cn_conv(i_cn)     = converged
        call flex_skew_kurt(p, np, cn_skews(i_cn), cn_kurts(i_cn))

        deallocate(p, p0)
    end do

    call rank_desc(cn_logls, cn_rank_logl)
    call rank_asc(cn_aics,   cn_rank_aic)
    call rank_asc(cn_bics,   cn_rank_bic)

    do i_cn = 1, cn_nmod
        idist = garch_dists(i_cn)
        if (idist == dist_normal) then
            write(*, cn_n) trim(cn_mnames(i_cn)), &
                cn_omegas(i_cn), cn_alphas(i_cn), cn_betas(i_cn), &
                cn_rhos(i_cn), cn_phis(i_cn), cn_thetas(i_cn), &
                cn_vol_anns(i_cn), cn_logls(i_cn), cn_aics(i_cn), cn_bics(i_cn), &
                cn_niters(i_cn), cn_skews(i_cn), cn_kurts(i_cn), &
                cn_rank_logl(i_cn), cn_rank_aic(i_cn), cn_rank_bic(i_cn)
        else if (idist == dist_t) then
            write(*, cn_t) trim(cn_mnames(i_cn)), &
                cn_omegas(i_cn), cn_alphas(i_cn), cn_betas(i_cn), &
                cn_rhos(i_cn), cn_phis(i_cn), cn_thetas(i_cn), cn_nus(i_cn), &
                cn_vol_anns(i_cn), cn_logls(i_cn), cn_aics(i_cn), cn_bics(i_cn), &
                cn_niters(i_cn), cn_skews(i_cn), cn_kurts(i_cn), &
                cn_rank_logl(i_cn), cn_rank_aic(i_cn), cn_rank_bic(i_cn)
        else
            write(*, cn_sk) trim(cn_mnames(i_cn)), &
                cn_omegas(i_cn), cn_alphas(i_cn), cn_betas(i_cn), &
                cn_rhos(i_cn), cn_phis(i_cn), cn_thetas(i_cn), &
                cn_nus(i_cn), cn_gammas(i_cn), &
                cn_vol_anns(i_cn), cn_logls(i_cn), cn_aics(i_cn), cn_bics(i_cn), &
                cn_niters(i_cn), cn_skews(i_cn), cn_kurts(i_cn), &
                cn_rank_logl(i_cn), cn_rank_aic(i_cn), cn_rank_bic(i_cn)
        end if
        if (.not. cn_conv(i_cn)) write(*, '(4X,A)') "WARNING: did not converge"
    end do

    ! ── Family GARCH models ──────────────────────────────────────────────────
    if (fit_fgarch) then
    write(*, '(/,A)') "  Family GARCH models (Hentschel 1995):"
    write(*, fg_hdr) "Model", "omega", "alpha", "beta", "lam", "nu", "b", "c", &
        "nu_t", "vol_ann%", "logL", "AIC", "BIC", "iter", "skew", "ekurt", "rklogL", "rkAIC", "rkBIC"
    write(*, '(A)') "  " // repeat("-", 155)

    call fgarch_set_data(ret, nobs)

    do i_fg = 1, fg_nmod
        call fgarch_set_dist(fg_dists(i_fg))
        np = fgarch_np()
        allocate(p(np), p0(np))

        ! Starting values: AGARCH point (lam=nu=1), slight negative shift for asymmetry.
        ! omega chosen so unconditional E[sigma] ~ s0 for lam=nu=1, b=0, c=0.
        call fgarch_inv_transform( &
            ret_std * (1.0_dp - 0.85_dp - 0.10_dp * 0.7979_dp), &
            0.10_dp, 0.85_dp, 1.0_dp, 1.0_dp, -0.1_dp, 0.0_dp, &
            8.0_dp, p0_fg)
        p0 = p0_fg(1:np)

        p = p0
        call bfgs_minimize(fgarch_obj, p, np, max_iter, gtol, fopt, niter, converged)

        call fgarch_transform(p, np, fg_omega, fg_alpha, fg_beta, fg_lam, fg_nu, fg_b, fg_c, fg_nu_t)
        call fgarch_vol_ann(p, np, vol_ann)
        logl = -nobs * fopt
        aic  = 2.0_dp*np             - 2.0_dp*logl
        bic  = np*log(real(nobs,dp)) - 2.0_dp*logl

        fg_omegas(i_fg)   = fg_omega
        fg_alphas(i_fg)   = fg_alpha
        fg_betas(i_fg)    = fg_beta
        fg_lams(i_fg)     = fg_lam
        fg_nus(i_fg)      = fg_nu
        fg_bs(i_fg)       = fg_b
        fg_cs(i_fg)       = fg_c
        fg_nu_ts(i_fg)    = fg_nu_t
        fg_vol_anns(i_fg) = vol_ann
        fg_logls(i_fg)    = logl
        fg_aics(i_fg)     = aic
        fg_bics(i_fg)     = bic
        fg_nps(i_fg)      = np
        fg_niters(i_fg)   = niter
        fg_conv(i_fg)     = converged
        call fgarch_skew_kurt(p, np, fg_skews(i_fg), fg_kurts(i_fg))

        deallocate(p, p0)
    end do

    call rank_desc(fg_logls, fg_rank_logl)
    call rank_asc(fg_aics,   fg_rank_aic)
    call rank_asc(fg_bics,   fg_rank_bic)

    do i_fg = 1, fg_nmod
        if (fg_dists(i_fg) == fg_dist_normal) then
            write(*, fg_fn) trim(fg_mnames(i_fg)), &
                fg_omegas(i_fg), fg_alphas(i_fg), fg_betas(i_fg), &
                fg_lams(i_fg), fg_nus(i_fg), fg_bs(i_fg), fg_cs(i_fg), &
                fg_vol_anns(i_fg), fg_logls(i_fg), fg_aics(i_fg), fg_bics(i_fg), &
                fg_niters(i_fg), fg_skews(i_fg), fg_kurts(i_fg), &
                fg_rank_logl(i_fg), fg_rank_aic(i_fg), fg_rank_bic(i_fg)
        else
            write(*, fg_ft) trim(fg_mnames(i_fg)), &
                fg_omegas(i_fg), fg_alphas(i_fg), fg_betas(i_fg), &
                fg_lams(i_fg), fg_nus(i_fg), fg_bs(i_fg), fg_cs(i_fg), fg_nu_ts(i_fg), &
                fg_vol_anns(i_fg), fg_logls(i_fg), fg_aics(i_fg), fg_bics(i_fg), &
                fg_niters(i_fg), fg_skews(i_fg), fg_kurts(i_fg), &
                fg_rank_logl(i_fg), fg_rank_aic(i_fg), fg_rank_bic(i_fg)
        end if
        if (.not. fg_conv(i_fg)) write(*, '(4X,A)') "WARNING: did not converge"
    end do
    end if   ! fit_fgarch

    ! ── Combined ranking ─────────────────────────────────────────────────────
    write(*, '(/,A)') "  Combined ranking by predictive logL:"
    write(*, comb_hdr) "Rank", "Model", "logL*", "AIC*", "BIC*", "skew", "ekurt"
    write(*, '(A)') "  " // repeat("-", 70)

    do iproc = proc_sv, proc_sv_lev
        do idist = dist_normal, dist_t
            imod = (iproc-1)*sv_n_dist + idist
            all_names(imod)  = sv_model_names(iproc, idist)
            all_logls(imod)  = sv_pred_logls(imod)
            all_aics(imod)   = 2.0_dp*sv_nps(imod)             - 2.0_dp*sv_pred_logls(imod)
            all_bics(imod)   = sv_nps(imod)*log(real(nobs,dp)) - 2.0_dp*sv_pred_logls(imod)
            all_skews(imod)  = sv_skews(imod)
            all_kurts(imod)  = sv_kurts(imod)
        end do
    end do
    do i_gp = 1, n_garch_procs
        do i_gd = 1, n_garch_dist
            imod = (i_gp - 1)*n_garch_dist + i_gd
            all_names(sv_nmod + imod)  = garch_mnames(i_gp, i_gd)
            all_logls(sv_nmod + imod)  = g_logls(imod)
            all_aics(sv_nmod + imod)   = g_aics(imod)
            all_bics(sv_nmod + imod)   = g_bics(imod)
            all_skews(sv_nmod + imod)  = g_skews(imod)
            all_kurts(sv_nmod + imod)  = g_kurts(imod)
        end do
    end do

    do i_cn = 1, cn_nmod
        imod = sv_nmod + garch_nmod + i_cn
        all_names(imod)  = cn_mnames(i_cn)
        all_logls(imod)  = cn_logls(i_cn)
        all_aics(imod)   = cn_aics(i_cn)
        all_bics(imod)   = cn_bics(i_cn)
        all_skews(imod)  = cn_skews(i_cn)
        all_kurts(imod)  = cn_kurts(i_cn)
    end do

    if (fit_fgarch) then
        do i_fg = 1, fg_nmod
            imod = sv_nmod + garch_nmod + cn_nmod + i_fg
            all_names(imod)  = fg_mnames(i_fg)
            all_logls(imod)  = fg_logls(i_fg)
            all_aics(imod)   = fg_aics(i_fg)
            all_bics(imod)   = fg_bics(i_fg)
            all_skews(imod)  = fg_skews(i_fg)
            all_kurts(imod)  = fg_kurts(i_fg)
        end do
    end if

    call rank_desc(all_logls, all_ranks)

    do irank = 1, nmod_total
        do j = 1, nmod_total
            if (all_ranks(j) == irank) then
                write(*, comb_fmt) irank, trim(all_names(j)), &
                    all_logls(j), all_aics(j), all_bics(j), all_skews(j), all_kurts(j)
                exit
            end if
        end do
    end do

    write(*,*)

end do

write(*, '(A)') "SV vol_ann% = sqrt(252 * exp(mu + sig_eta^2 / (2*(1-phi^2)))) * 100"
write(*, '(A)') "par3: theta (NAGARCH), gamma (GJR/EGARCH)"
write(*, '(A)') "logL* for SV: sum_t log p(y_t | Y_{t-1}) using KF predicted state (exp(a_t + P_t/2) as variance)"
write(*, '(A)') "logL* for GARCH: standard GARCH log-likelihood (identical scale -- directly comparable)"
write(*, '(A)') "logL* for FGARCH: same scale as GARCH (filter uses past sigma to predict current y)"
write(*, '(A)') "FGARCH: h_t=(sigma_t^lam-1)/lam = omega + alpha*sigma_{t-1}^nu*f(z)^nu + beta*h_{t-1}"
write(*, '(A)') "        f(z) = |z-b| - c*(z-b),  z = y/sigma  (b=shift, c=rotation of news impact curve)"
write(*, '(A)') "Abs-value models (AVGARCH/ANAGARCH/AGJR): update sigma_t (not variance h_t)"
write(*, '(A)') "  AVGARCH:  sigma_t = omega + alpha*|y_{t-1}| + beta*sigma_{t-1}"
write(*, '(A)') "  ANAGARCH: sigma_t = omega + alpha*|y_{t-1}-theta*sigma_{t-1}| + beta*sigma_{t-1}"
write(*, '(A)') "  AGJR:     sigma_t = omega + (alpha+gamma*I_{t-1})*|y_{t-1}| + beta*sigma_{t-1}"
write(*, '(A)') "  vol_ann% = omega/(1-alpha-beta) * sqrt(252) * 100  (approx unconditional sigma)"
write(*, '(A)') "Component NAGARCH (Engle-Lee 1999):"
write(*, '(A)') "  q_t = omega + rho*q_{t-1} + phi*(r²_{t-1} - h_{t-1})"
write(*, '(A)') "  h_t = q_t + alpha*(r²_{t-1} - q_{t-1}) + beta*(h_{t-1} - q_{t-1})"
write(*, '(A)') "  r_t = y_t - theta*sqrt(h_t)  (NAGARCH shifted residual)"
write(*, '(A)') "  vol_ann% = sqrt(252 * h_bar) * 100,  h_bar = q_bar*(1-alpha-beta)/D1"
write(*, '(A)') "  q_bar = omega/((1-rho) - phi*theta^2*(1-alpha-beta)/D1),  D1 = 1-alpha*(1+theta^2)-beta"

call cpu_time(t_end)
write(*, '(/,A,F8.3,A)') "elapsed time: ", t_end - t_start, " s"

end program xfit_sv_garch_returns

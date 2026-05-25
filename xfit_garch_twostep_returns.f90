! Two-step GARCH/distribution fitting for demeaned log returns.
!
! Usage:
!   xfit_garch_twostep_returns.exe [prices_file] [output_csv]
!
! This is a fast screening approximation to joint GARCH/distribution maximum
! likelihood. For each asset and GARCH model, the program first fits the
! volatility model once under normal innovations. It then computes standardized
! residuals z_t = r_t / sqrt(h_t), fits each supported iid innovation
! distribution to those fixed residuals, and scores the pair using the fixed
! volatility path:
!
!   log L = sum_t log f(z_t) - 0.5 * sum_t log(h_t)
!
! Because the volatility parameters are not re-optimized for each distribution,
! these likelihoods are generally lower than the corresponding joint fits from
! xfit_gen_garch_dist_returns.f90. Use this program for fast distribution
! screening and the joint-fit driver for final likelihood comparisons.
!
module garch_twostep_returns_mod
    use kind_mod, only: dp
    use csv_mod, only: read_price_csv, print_price_sample_info
    use stats_mod, only: mean, variance
    use distributions_mod, only: pdf_normal, pdf_t, pdf_ged, pdf_logistic, pdf_laplace, pdf_sech, &
        pdf_nig, pdf_fs_skewt
    use garch_fit_dist_mod, only: garch_dist_fit_result_t, fit_garch_dist_model, garch_dist_variance_path, &
        garch_dist_persist, model_param_count, dist_param_count, dist_nu_value, dist_xi_value, dist_alpha_value
    use bfgs_mod, only: bfgs_minimize
    implicit none
    private

    public :: run_garch_twostep_returns

    integer, parameter :: max_assets = 10**6
    integer, parameter :: max_iter = 100
    real(dp), parameter :: gtol = 1.0e-7_dp
    real(dp), parameter :: min_h = 1.0e-12_dp
    real(dp), parameter :: min_pdf = 1.0e-300_dp
    real(dp), parameter :: trading_days = 252.0_dp
    character(len=*), parameter :: default_prices_file = "spy_efa_eem_tlt_lqd.csv"
    character(len=16), parameter :: models(*) = [character(len=16) :: "SYMM_GARCH", "NAGARCH"]
    character(len=10), parameter :: dists(*) = [character(len=10) :: &
        "NORMAL", "T"] ! , "SECH", "GED", "LAPLACE", "LOGISTIC", "NIG", "FS_SKEWT"]
    integer, parameter :: n_model = size(models)
    integer, parameter :: n_dist = size(dists)
    integer, parameter :: n_combo = n_model*n_dist

    real(dp), allocatable, save :: obj_z(:)
    character(len=10), save :: obj_dist = ""
    integer, save :: obj_np = 0

    type :: twostep_row_t
        character(len=16) :: model = ""
        character(len=10) :: dist = ""
        character(len=32) :: asset = ""
        type(garch_dist_fit_result_t) :: garch_fit
        real(dp) :: shape = 0.0_dp
        real(dp) :: shape2 = 0.0_dp
        real(dp) :: persist = 0.0_dp
        real(dp) :: vol_ann = 0.0_dp
        real(dp) :: logl = -huge(1.0_dp)
        real(dp) :: aic = huge(1.0_dp)
        real(dp) :: bic = huge(1.0_dp)
        real(dp) :: fit_sec = 0.0_dp
        integer :: niter = 0
        logical :: converged = .false.
    end type twostep_row_t

contains

    subroutine run_garch_twostep_returns()
        integer, allocatable :: dates(:)
        character(len=32), allocatable :: col_names(:)
        real(dp), allocatable :: prices(:,:), ret(:), h(:), z(:), asset_fit_seconds(:)
        integer, allocatable :: asset_converged_count(:), asset_iter_total(:)
        type(twostep_row_t) :: rows(n_combo)
        type(garch_dist_fit_result_t) :: normal_fit
        integer :: row_aic_rank(n_combo), row_bic_rank(n_combo)
        integer :: row_model_idx(n_combo), row_dist_idx(n_combo)
        integer :: aic_wins(n_model,n_dist), bic_wins(n_model,n_dist)
        integer :: aic_rank_sum(n_model,n_dist), bic_rank_sum(n_model,n_dist), rank_count(n_model,n_dist)
        real(dp) :: fit_seconds(n_model,n_dist)
        integer :: nprices, ncols, nobs, icol, imod, idist, ifit, nfit
        integer :: nparam, clock_start, clock_end, clock_rate, fit_clock_start, fit_clock_end
        integer :: csv_unit
        real(dp) :: ret_mean, elapsed_s, vol_ann, dist_logl, h_log_sum, garch_fit_sec, dist_only_sec
        character(len=256) :: prices_file, output_csv
        logical :: write_csv
        print*,"default prices file: " // trim(default_prices_file)
        prices_file = default_prices_file
        output_csv = ""
        if (command_argument_count() >= 1) call get_command_argument(1, prices_file)
        if (command_argument_count() >= 2) call get_command_argument(2, output_csv)
        write_csv = len_trim(output_csv) > 0
        print*,"write_csv =",write_csv ! debug
        call system_clock(clock_start, clock_rate)
        aic_wins = 0
        bic_wins = 0
        aic_rank_sum = 0
        bic_rank_sum = 0
        rank_count = 0
        fit_seconds = 0.0_dp

!        stop "before calling read_price_csv" ! debug
        call read_price_csv(prices_file, dates, col_names, prices, max_col=max_assets)
        nprices = size(prices, 1)
        ncols = size(prices, 2)
        nobs = nprices - 1
        print*,"ncols, nobs =", ncols, nobs ! debug
        allocate(ret(nobs), h(nobs), z(nobs))
        allocate(asset_fit_seconds(ncols), asset_converged_count(ncols), asset_iter_total(ncols))
        asset_fit_seconds = 0.0_dp
        asset_converged_count = 0
        asset_iter_total = 0

        call print_price_sample_info(trim(prices_file), dates, ncols)
        write(*,'(A,I0)') "Maximum assets read: ", max_assets
        if (write_csv) then
            open(newunit=csv_unit, file=trim(output_csv), status="replace", action="write")
            write(csv_unit,'(A)') "method,model,dist,asset,omega,alpha,gamma,beta,theta,nu,xi,dist_alpha,persist,vol_ann_pct,logL,AIC,BIC,iter,conv,AIC_rank,BIC_rank,fit_sec"
        end if
        write(*,'(/,A)') "Two-step GARCH/distribution fits: normal GARCH volatility, then iid distribution fit to standardized residuals"
        write(*,'(A)') "Method       Model            Dist       Asset        omega   alpha   gamma    beta   theta      nu      xi  dist_alpha  persist  vol_ann%        logL         AIC         BIC  iter conv AIC_rank BIC_rank    fit_sec"
        write(*,'(A)') repeat("-", 210)

        do icol = 1, ncols
            ret = log(prices(2:nprices,icol) / prices(1:nprices-1,icol))
            ret_mean = mean(ret)
            ret = ret - ret_mean
            vol_ann = sqrt(max(variance(ret), 0.0_dp)*trading_days)*100.0_dp
            nfit = 0

            do imod = 1, n_model
                call system_clock(fit_clock_start)
                call fit_garch_dist_model(trim(models(imod)), "NORMAL", ret, max_iter, gtol, normal_fit)
                call garch_dist_variance_path(trim(models(imod)), ret, normal_fit%params, h)
                h = max(h, min_h)
                z = ret / sqrt(h)
                h_log_sum = sum(log(h))
                call system_clock(fit_clock_end)
                garch_fit_sec = real(fit_clock_end - fit_clock_start, dp) / real(clock_rate, dp)
                asset_fit_seconds(icol) = asset_fit_seconds(icol) + garch_fit_sec

                do idist = 1, n_dist
                    nfit = nfit + 1
                    row_model_idx(nfit) = imod
                    row_dist_idx(nfit) = idist
                    rows(nfit)%model = trim(models(imod))
                    rows(nfit)%dist = trim(dists(idist))
                    rows(nfit)%asset = trim(col_names(icol))
                    rows(nfit)%garch_fit = normal_fit
                    rows(nfit)%persist = garch_dist_persist(trim(models(imod)), normal_fit%params)
                    rows(nfit)%vol_ann = vol_ann

                    if (idist == 1) then
                        rows(nfit)%shape = 0.0_dp
                        rows(nfit)%shape2 = 0.0_dp
                        rows(nfit)%niter = normal_fit%niter
                        rows(nfit)%converged = normal_fit%converged
                        rows(nfit)%fit_sec = garch_fit_sec
                        dist_logl = standardized_loglik(trim(dists(idist)), z, rows(nfit)%shape, rows(nfit)%shape2)
                    else
                        call system_clock(fit_clock_start)
                        call fit_standardized_dist(trim(dists(idist)), z, rows(nfit)%shape, rows(nfit)%shape2, &
                                                   dist_logl, rows(nfit)%niter, rows(nfit)%converged)
                        call system_clock(fit_clock_end)
                        dist_only_sec = real(fit_clock_end - fit_clock_start, dp) / real(clock_rate, dp)
                        rows(nfit)%fit_sec = dist_only_sec
                        asset_fit_seconds(icol) = asset_fit_seconds(icol) + dist_only_sec
                    end if
                    rows(nfit)%logl = dist_logl - 0.5_dp*h_log_sum
                    nparam = model_param_count(trim(models(imod))) + dist_param_count(trim(dists(idist)))
                    rows(nfit)%aic = 2.0_dp*real(nparam, dp) - 2.0_dp*rows(nfit)%logl
                    rows(nfit)%bic = log(real(nobs, dp))*real(nparam, dp) - 2.0_dp*rows(nfit)%logl
                    fit_seconds(imod,idist) = fit_seconds(imod,idist) + rows(nfit)%fit_sec
                    if (rows(nfit)%converged) asset_converged_count(icol) = asset_converged_count(icol) + 1
                    asset_iter_total(icol) = asset_iter_total(icol) + rows(nfit)%niter
                end do
            end do

            call rank_rows(rows, nfit, row_aic_rank, row_bic_rank)
            do ifit = 1, nfit
                imod = row_model_idx(ifit)
                idist = row_dist_idx(ifit)
                if (row_aic_rank(ifit) == 1) aic_wins(imod,idist) = aic_wins(imod,idist) + 1
                if (row_bic_rank(ifit) == 1) bic_wins(imod,idist) = bic_wins(imod,idist) + 1
                aic_rank_sum(imod,idist) = aic_rank_sum(imod,idist) + row_aic_rank(ifit)
                bic_rank_sum(imod,idist) = bic_rank_sum(imod,idist) + row_bic_rank(ifit)
                rank_count(imod,idist) = rank_count(imod,idist) + 1
            end do

            do ifit = 1, nfit
                call print_row(rows(ifit), row_aic_rank(ifit), row_bic_rank(ifit), write_csv, csv_unit)
            end do
        end do

        if (write_csv) then
            close(csv_unit)
            write(*,'(/,A,A)') "Wrote fit table CSV: ", trim(output_csv)
        end if
        call print_asset_summary(col_names, asset_fit_seconds, asset_converged_count, asset_iter_total)
        call print_combo_summary(aic_wins, bic_wins, aic_rank_sum, bic_rank_sum, rank_count, fit_seconds)
        call system_clock(clock_end)
        elapsed_s = real(clock_end - clock_start, dp) / real(clock_rate, dp)
        write(*,'(/,A,F10.3,A)') "Elapsed wall time: ", elapsed_s, " seconds"
    end subroutine run_garch_twostep_returns

    subroutine fit_standardized_dist(dist_name, z, shape, shape2, logl, niter, converged)
        character(len=*), intent(in) :: dist_name
        real(dp), intent(in) :: z(:)
        real(dp), intent(out) :: shape, shape2, logl
        integer, intent(out) :: niter
        logical, intent(out) :: converged
        real(dp), allocatable :: p(:), grad(:)
        real(dp) :: f_best

        obj_dist = trim(dist_name)
        obj_np = dist_param_count(obj_dist)
        if (allocated(obj_z)) deallocate(obj_z)
        allocate(obj_z(size(z)))
        obj_z = z

        if (obj_np == 0) then
            shape = 0.0_dp
            shape2 = 0.0_dp
            logl = standardized_loglik(obj_dist, z, shape, shape2)
            niter = 0
            converged = .true.
            return
        end if

        allocate(p(obj_np), grad(obj_np))
        call initial_dist_params(obj_dist, p)
        call bfgs_minimize(dist_obj, p, obj_np, 500, gtol, f_best, niter, converged)
        call unpack_dist_params(obj_dist, p, shape, shape2)
        logl = -f_best*real(size(z), dp)
        deallocate(p, grad)
    end subroutine fit_standardized_dist

    subroutine dist_obj(p, np, f, g)
        integer, intent(in) :: np
        real(dp), intent(in) :: p(np)
        real(dp), intent(out) :: f, g(np)
        real(dp), allocatable :: pp(:), pm(:)
        real(dp) :: fp, fm, step
        integer :: j

        f = dist_nll_from_p(p, np)
        allocate(pp(np), pm(np))
        do j = 1, np
            step = 1.0e-5_dp*max(1.0_dp, abs(p(j)))
            pp = p
            pm = p
            pp(j) = pp(j) + step
            pm(j) = pm(j) - step
            fp = dist_nll_from_p(pp, np)
            fm = dist_nll_from_p(pm, np)
            g(j) = (fp - fm) / (2.0_dp*step)
        end do
        deallocate(pp, pm)
    end subroutine dist_obj

    real(dp) function dist_nll_from_p(p, np)
        integer, intent(in) :: np
        real(dp), intent(in) :: p(np)
        real(dp) :: shape, shape2

        call unpack_dist_params(obj_dist, p, shape, shape2)
        dist_nll_from_p = -standardized_loglik(obj_dist, obj_z, shape, shape2) / real(size(obj_z), dp)
        if (dist_nll_from_p /= dist_nll_from_p .or. dist_nll_from_p > 1.0e29_dp) dist_nll_from_p = 1.0e30_dp
    end function dist_nll_from_p

    subroutine initial_dist_params(dist_name, p)
        character(len=*), intent(in) :: dist_name
        real(dp), intent(out) :: p(:)

        select case (trim(dist_name))
        case ("T")
            p(1) = pack_bounded(8.0_dp, 2.01_dp, 100.0_dp)
        case ("GED")
            p(1) = log(1.5_dp)
        case ("NIG")
            p(1) = pack_bounded(3.0_dp, 0.1_dp, 20.0_dp)
        case ("FS_SKEWT")
            p(1) = pack_bounded(8.0_dp, 2.01_dp, 100.0_dp)
            p(2) = log(1.0_dp)
        case default
            if (size(p) > 0) p = 0.0_dp
        end select
    end subroutine initial_dist_params

    subroutine unpack_dist_params(dist_name, p, shape, shape2)
        character(len=*), intent(in) :: dist_name
        real(dp), intent(in) :: p(:)
        real(dp), intent(out) :: shape, shape2

        shape = 0.0_dp
        shape2 = 0.0_dp
        select case (trim(dist_name))
        case ("T")
            shape = unpack_bounded(p(1), 2.01_dp, 100.0_dp)
        case ("GED")
            shape = exp(clamp(p(1), -40.0_dp, 40.0_dp))
        case ("NIG")
            shape = unpack_bounded(p(1), 0.1_dp, 20.0_dp)
        case ("FS_SKEWT")
            shape = unpack_bounded(p(1), 2.01_dp, 100.0_dp)
            shape2 = exp(clamp(p(2), -40.0_dp, 40.0_dp))
        end select
    end subroutine unpack_dist_params

    real(dp) function standardized_loglik(dist_name, z, shape, shape2)
        character(len=*), intent(in) :: dist_name
        real(dp), intent(in) :: z(:), shape, shape2
        integer :: i

        standardized_loglik = 0.0_dp
        do i = 1, size(z)
            standardized_loglik = standardized_loglik + log(max(dist_pdf(trim(dist_name), z(i), shape, shape2), min_pdf))
        end do
    end function standardized_loglik

    real(dp) function dist_pdf(dist_name, z, shape, shape2)
        character(len=*), intent(in) :: dist_name
        real(dp), intent(in) :: z, shape, shape2

        select case (trim(dist_name))
        case ("NORMAL")
            dist_pdf = pdf_normal(z)
        case ("T")
            dist_pdf = pdf_t(z, shape)
        case ("SECH")
            dist_pdf = pdf_sech(z)
        case ("GED")
            dist_pdf = pdf_ged(z, shape)
        case ("LAPLACE")
            dist_pdf = pdf_laplace(z)
        case ("LOGISTIC")
            dist_pdf = pdf_logistic(z)
        case ("NIG")
            dist_pdf = pdf_nig(z, shape)
        case ("FS_SKEWT")
            dist_pdf = pdf_fs_skewt(z, shape, shape2)
        case default
            dist_pdf = min_pdf
        end select
    end function dist_pdf

    subroutine rank_rows(rows, nfit, aic_rank, bic_rank)
        type(twostep_row_t), intent(in) :: rows(:)
        integer, intent(in) :: nfit
        integer, intent(out) :: aic_rank(:), bic_rank(:)
        integer :: i, j

        do i = 1, nfit
            aic_rank(i) = 1
            bic_rank(i) = 1
            do j = 1, nfit
                if (rows(j)%aic < rows(i)%aic) aic_rank(i) = aic_rank(i) + 1
                if (rows(j)%bic < rows(i)%bic) bic_rank(i) = bic_rank(i) + 1
            end do
        end do
    end subroutine rank_rows

    subroutine print_row(row, aic_rank, bic_rank, write_csv, csv_unit)
        type(twostep_row_t), intent(in) :: row
        integer, intent(in) :: aic_rank, bic_rank, csv_unit
        logical, intent(in) :: write_csv

        write(*,'(A10,1X,A16,1X,A10,1X,A9,ES12.3,4F8.4,1X,A8,1X,A8,1X,A10,F9.4,F10.2,F12.2,F12.2,F12.2,I6,1X,L1,2I9,1X,F10.3)') &
            "TWO_STEP", trim(row%model), trim(row%dist), trim(row%asset), row%garch_fit%params%omega, &
            row%garch_fit%params%alpha, row%garch_fit%params%gamma, row%garch_fit%params%beta, &
            row%garch_fit%params%theta, dist_nu_value(trim(row%dist), row%shape), &
            dist_xi_value(trim(row%dist), row%shape2), dist_alpha_value(trim(row%dist), row%shape), &
            row%persist, row%vol_ann, row%logl, row%aic, row%bic, row%niter, row%converged, &
            aic_rank, bic_rank, row%fit_sec
        if (write_csv) then
            write(csv_unit,'(A,",",A,",",A,",",A,",",ES16.8,",",F12.6,",",F12.6,",",F12.6,",",F12.6,",",A,",",A,",",A,",",F12.6,",",F12.6,",",F16.6,",",F16.6,",",F16.6,",",I0,",",L1,",",I0,",",I0,",",F12.6)') &
                "TWO_STEP", trim(row%model), trim(row%dist), trim(row%asset), row%garch_fit%params%omega, &
                row%garch_fit%params%alpha, row%garch_fit%params%gamma, row%garch_fit%params%beta, &
                row%garch_fit%params%theta, dist_nu_value(trim(row%dist), row%shape), &
                dist_xi_value(trim(row%dist), row%shape2), dist_alpha_value(trim(row%dist), row%shape), &
                row%persist, row%vol_ann, row%logl, row%aic, row%bic, row%niter, row%converged, &
                aic_rank, bic_rank, row%fit_sec
        end if
    end subroutine print_row

    subroutine print_asset_summary(asset_names, asset_seconds, converged_count, iter_total)
        character(len=*), intent(in) :: asset_names(:)
        real(dp), intent(in) :: asset_seconds(:)
        integer, intent(in) :: converged_count(:), iter_total(:)
        integer :: ia

        write(*,'(/,A)') "Fit time by asset"
        write(*,'(A)') repeat("-", 58)
        write(*,'(A12,1X,A12,1X,A11,1X,A10)') "Asset", "fit_sec", "converged", "iter_total"
        write(*,'(A)') repeat("-", 58)
        do ia = 1, size(asset_names)
            write(*,'(A12,1X,F12.3,1X,I11,1X,I10)') trim(asset_names(ia)), asset_seconds(ia), &
                converged_count(ia), iter_total(ia)
        end do
        write(*,'(A)') repeat("-", 58)
    end subroutine print_asset_summary

    subroutine print_combo_summary(aic_wins, bic_wins, aic_rank_sum, bic_rank_sum, rank_count, fit_seconds)
        integer, intent(in) :: aic_wins(n_model,n_dist), bic_wins(n_model,n_dist)
        integer, intent(in) :: aic_rank_sum(n_model,n_dist), bic_rank_sum(n_model,n_dist)
        integer, intent(in) :: rank_count(n_model,n_dist)
        real(dp), intent(in) :: fit_seconds(n_model,n_dist)
        real(dp) :: mean_aic_rank, mean_bic_rank, freq
        integer :: im, id, model_total, dist_total, grand_total

        write(*,'(/,A)') "Two-step model-distribution selection summary"
        write(*,'(A)') repeat("-", 86)
        write(*,'(A16,1X,A10,1X,A8,1X,A8,1X,A13,1X,A13,1X,A10)') &
            "Model", "Dist", "AIC_wins", "BIC_wins", "mean_AIC_rank", "mean_BIC_rank", "fit_sec"
        write(*,'(A)') repeat("-", 86)
        do im = 1, n_model
            do id = 1, n_dist
                if (rank_count(im,id) > 0) then
                    mean_aic_rank = real(aic_rank_sum(im,id), dp) / real(rank_count(im,id), dp)
                    mean_bic_rank = real(bic_rank_sum(im,id), dp) / real(rank_count(im,id), dp)
                else
                    mean_aic_rank = 0.0_dp
                    mean_bic_rank = 0.0_dp
                end if
                write(*,'(A16,1X,A10,1X,I8,1X,I8,1X,F13.3,1X,F13.3,1X,F10.3)') &
                    trim(models(im)), trim(dists(id)), aic_wins(im,id), bic_wins(im,id), &
                    mean_aic_rank, mean_bic_rank, fit_seconds(im,id)
            end do
        end do
        write(*,'(A)') repeat("-", 86)

        write(*,'(/,A)') "AIC-selected distribution frequencies within each GARCH model"
        write(*,'(A)') repeat("-", 58)
        write(*,'(A16,1X,A10,1X,A9,1X,A9)') "Model", "Dist", "AIC_count", "frequency"
        write(*,'(A)') repeat("-", 58)
        do im = 1, n_model
            model_total = sum(aic_wins(im,:))
            do id = 1, n_dist
                freq = 0.0_dp
                if (model_total > 0) freq = real(aic_wins(im,id), dp) / real(model_total, dp)
                write(*,'(A16,1X,A10,1X,I9,1X,F9.3)') trim(models(im)), trim(dists(id)), aic_wins(im,id), freq
            end do
        end do
        write(*,'(A)') repeat("-", 58)

        grand_total = sum(aic_wins)
        write(*,'(/,A)') "AIC-selected distribution frequencies across GARCH models"
        write(*,'(A)') repeat("-", 42)
        write(*,'(A10,1X,A9,1X,A9)') "Dist", "AIC_count", "frequency"
        write(*,'(A)') repeat("-", 42)
        do id = 1, n_dist
            dist_total = sum(aic_wins(:,id))
            freq = 0.0_dp
            if (grand_total > 0) freq = real(dist_total, dp) / real(grand_total, dp)
            write(*,'(A10,1X,I9,1X,F9.3)') trim(dists(id)), dist_total, freq
        end do
        write(*,'(A)') repeat("-", 42)
    end subroutine print_combo_summary

    real(dp) function pack_bounded(x, lo, hi)
        real(dp), intent(in) :: x, lo, hi
        real(dp) :: xx

        xx = min(max(x, lo + 1.0e-8_dp), hi - 1.0e-8_dp)
        pack_bounded = log((xx - lo) / (hi - xx))
    end function pack_bounded

    real(dp) function unpack_bounded(q, lo, hi)
        real(dp), intent(in) :: q, lo, hi

        unpack_bounded = lo + (hi - lo) / (1.0_dp + exp(-clamp(q, -40.0_dp, 40.0_dp)))
    end function unpack_bounded

    real(dp) function clamp(x, lo, hi)
        real(dp), intent(in) :: x, lo, hi

        clamp = min(max(x, lo), hi)
    end function clamp

end module garch_twostep_returns_mod

program xfit_garch_twostep_returns
    use garch_twostep_returns_mod, only: run_garch_twostep_returns
    implicit none

    call run_garch_twostep_returns()
end program xfit_garch_twostep_returns

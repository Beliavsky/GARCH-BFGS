! Fit selected GARCH-family models with selected innovation distributions
! to demeaned log returns computed from a price CSV.
!
! Usage:
!   xfit_gen_garch_dist_returns.exe [prices_file] [output_csv]

program xfit_gen_garch_dist_returns
    use kind_mod, only: dp
    use csv_mod, only: read_price_csv, print_price_sample_info
    use stats_mod, only: mean, sd, variance
    use garch_fit_dist_mod, only: garch_dist_fit_result_t, fit_garch_dist_model, &
        model_param_count, dist_param_count, dist_nu_value, dist_xi_value, dist_alpha_value
    implicit none

    character(len=*), parameter :: default_prices_file = "spy_efa_eem_tlt_lqd.csv"
    real(dp), parameter :: trading_days = 252.0_dp
    integer, parameter :: max_assets = 3, max_iter = 100
    real(dp), parameter :: gtol = 1.0e-7_dp
    character(len=16), parameter :: models(*) = [character(len=16) :: &
        "SYMM_GARCH", "NAGARCH"] ! , "GJR_GARCH", "GJR_SIGNED", "EGARCH", "QGARCH"]
    character(len=10), parameter :: dists(*) = [character(len=10) :: &
        "NORMAL", "T", "SECH", "GED", "LAPLACE", "LOGISTIC", "NIG", "FS_SKEWT"]
    integer, parameter :: n_model = size(models)
    integer, parameter :: n_dist = size(dists)
    integer, parameter :: n_combo = n_model*n_dist

    integer, allocatable :: dates(:)
    character(len=32), allocatable :: col_names(:)
    real(dp), allocatable :: prices(:,:), ret(:), asset_fit_seconds(:)
    type(garch_dist_fit_result_t) :: rows(n_combo)
    integer :: row_aic_rank(n_combo), row_bic_rank(n_combo)
    integer :: row_model_idx(n_combo), row_dist_idx(n_combo)
    integer, allocatable :: asset_converged_count(:), asset_iter_total(:)
    integer :: aic_wins(n_model,n_dist), bic_wins(n_model,n_dist)
    integer :: aic_rank_sum(n_model,n_dist), bic_rank_sum(n_model,n_dist), rank_count(n_model,n_dist)
    real(dp) :: fit_seconds(n_model,n_dist)
    integer :: nprices, ncols, nobs, icol, imod, idist, ifit, jfit, nfit
    integer :: nparam, clock_start, clock_end, clock_rate, fit_clock_start, fit_clock_end
    integer :: csv_unit
    real(dp) :: ret_mean, ret_std, elapsed_s, logl, aic, bic, vol_ann
    real(dp) :: row_aic(n_combo), row_bic(n_combo), row_logl(n_combo), row_vol_ann(n_combo), row_fit_sec(n_combo)
    character(len=256) :: prices_file, output_csv
    logical :: write_csv

    prices_file = default_prices_file
    output_csv = ""
    if (command_argument_count() >= 1) call get_command_argument(1, prices_file)
    if (command_argument_count() >= 2) call get_command_argument(2, output_csv)
    write_csv = len_trim(output_csv) > 0

    call system_clock(clock_start, clock_rate)
    aic_wins = 0
    bic_wins = 0
    aic_rank_sum = 0
    bic_rank_sum = 0
    rank_count = 0
    fit_seconds = 0.0_dp
    print*,"max_assets =",max_assets
    call read_price_csv(trim(prices_file), dates, col_names, prices, max_col=max_assets)
    nprices = size(prices, 1)
    ncols = size(prices, 2)
    nobs = nprices - 1
    allocate(ret(nobs))
    allocate(asset_fit_seconds(ncols))
    allocate(asset_converged_count(ncols), asset_iter_total(ncols))
    asset_fit_seconds = 0.0_dp
    asset_converged_count = 0
    asset_iter_total = 0

    call print_price_sample_info(trim(prices_file), dates, ncols)
    if (write_csv) then
        open(newunit=csv_unit, file=trim(output_csv), status="replace", action="write")
        write(csv_unit,'(A)') "model,dist,asset,omega,alpha,gamma,beta,theta,nu,xi,dist_alpha,persist,vol_ann_pct,logL,AIC,BIC,iter,conv,AIC_rank,BIC_rank,fit_sec"
    end if
    write(*,'(A)') "Model            Dist       Asset        omega   alpha   gamma    beta   theta      nu      xi  dist_alpha  persist  vol_ann%        logL         AIC         BIC  iter conv AIC_rank BIC_rank    fit_sec"
    write(*,'(A)') repeat("-", 197)
    do icol = 1, ncols
        ret = log(prices(2:nprices,icol) / prices(1:nprices-1,icol))
        ret_mean = mean(ret)
        ret = ret - ret_mean
        ret_std = sd(ret)
        nfit = 0

        do imod = 1, n_model
            do idist = 1, n_dist
                nfit = nfit + 1
                row_model_idx(nfit) = imod
                row_dist_idx(nfit) = idist

                call system_clock(fit_clock_start)
                call fit_garch_dist_model(trim(models(imod)), trim(dists(idist)), ret, max_iter, gtol, rows(nfit))
                call system_clock(fit_clock_end)
                row_fit_sec(nfit) = real(fit_clock_end - fit_clock_start, dp) / real(clock_rate, dp)
                fit_seconds(imod,idist) = fit_seconds(imod,idist) + &
                    row_fit_sec(nfit)
                asset_fit_seconds(icol) = asset_fit_seconds(icol) + row_fit_sec(nfit)
                if (rows(nfit)%converged) asset_converged_count(icol) = asset_converged_count(icol) + 1
                asset_iter_total(icol) = asset_iter_total(icol) + rows(nfit)%niter

                nparam = model_param_count(trim(models(imod))) + dist_param_count(trim(dists(idist)))
                logl = -real(nobs, dp)*rows(nfit)%nll
                aic = 2.0_dp*real(nparam, dp) - 2.0_dp*logl
                bic = log(real(nobs, dp))*real(nparam, dp) - 2.0_dp*logl
                vol_ann = sqrt(max(variance(ret), 0.0_dp)*trading_days)*100.0_dp

                row_logl(nfit) = logl
                row_aic(nfit) = aic
                row_bic(nfit) = bic
                row_vol_ann(nfit) = vol_ann
            end do
        end do

        do ifit = 1, nfit
            row_aic_rank(ifit) = 1
            row_bic_rank(ifit) = 1
            do jfit = 1, nfit
                if (row_aic(jfit) < row_aic(ifit)) row_aic_rank(ifit) = row_aic_rank(ifit) + 1
                if (row_bic(jfit) < row_bic(ifit)) row_bic_rank(ifit) = row_bic_rank(ifit) + 1
            end do
            imod = row_model_idx(ifit)
            idist = row_dist_idx(ifit)
            if (row_aic_rank(ifit) == 1) aic_wins(imod,idist) = aic_wins(imod,idist) + 1
            if (row_bic_rank(ifit) == 1) bic_wins(imod,idist) = bic_wins(imod,idist) + 1
            aic_rank_sum(imod,idist) = aic_rank_sum(imod,idist) + row_aic_rank(ifit)
            bic_rank_sum(imod,idist) = bic_rank_sum(imod,idist) + row_bic_rank(ifit)
            rank_count(imod,idist) = rank_count(imod,idist) + 1
        end do

        do ifit = 1, nfit
            write(*,'(A16,1X,A10,1X,A9,ES12.3,4F8.4,1X,A8,1X,A8,1X,A10,F9.4,F10.2,F12.2,F12.2,F12.2,I6,1X,L1,2I9,1X,F10.3)') &
                trim(rows(ifit)%model), trim(rows(ifit)%dist), trim(col_names(icol)), &
                rows(ifit)%params%omega, rows(ifit)%params%alpha, rows(ifit)%params%gamma, &
                rows(ifit)%params%beta, rows(ifit)%params%theta, dist_nu_value(trim(rows(ifit)%dist), rows(ifit)%shape), &
                dist_xi_value(trim(rows(ifit)%dist), rows(ifit)%shape2), &
                dist_alpha_value(trim(rows(ifit)%dist), rows(ifit)%shape), rows(ifit)%persist, row_vol_ann(ifit), &
                row_logl(ifit), row_aic(ifit), row_bic(ifit), rows(ifit)%niter, rows(ifit)%converged, &
                row_aic_rank(ifit), row_bic_rank(ifit), row_fit_sec(ifit)
            if (write_csv) then
                write(csv_unit,'(A,",",A,",",A,",",ES16.8,",",F12.6,",",F12.6,",",F12.6,",",F12.6,",",A,",",A,",",A,",",F12.6,",",F12.6,",",F16.6,",",F16.6,",",F16.6,",",I0,",",L1,",",I0,",",I0,",",F12.6)') &
                    trim(rows(ifit)%model), trim(rows(ifit)%dist), trim(col_names(icol)), &
                    rows(ifit)%params%omega, rows(ifit)%params%alpha, rows(ifit)%params%gamma, &
                    rows(ifit)%params%beta, rows(ifit)%params%theta, &
                    dist_nu_value(trim(rows(ifit)%dist), rows(ifit)%shape), &
                    dist_xi_value(trim(rows(ifit)%dist), rows(ifit)%shape2), &
                    dist_alpha_value(trim(rows(ifit)%dist), rows(ifit)%shape), rows(ifit)%persist, row_vol_ann(ifit), &
                    row_logl(ifit), row_aic(ifit), row_bic(ifit), rows(ifit)%niter, rows(ifit)%converged, &
                    row_aic_rank(ifit), row_bic_rank(ifit), row_fit_sec(ifit)
            end if
        end do
    end do
    if (write_csv) then
        close(csv_unit)
        write(*,'(/,A,A)') "Wrote fit table CSV: ", trim(output_csv)
    end if
    call print_asset_time_summary(col_names, asset_fit_seconds, asset_converged_count, asset_iter_total)
    call print_combo_summary(aic_wins, bic_wins, aic_rank_sum, bic_rank_sum, rank_count, fit_seconds)
    call system_clock(clock_end)
    elapsed_s = real(clock_end - clock_start, dp) / real(clock_rate, dp)
    write(*,'(/,A,F10.3,A)') "Elapsed wall time: ", elapsed_s, " seconds"

contains

    subroutine print_asset_time_summary(asset_names, asset_seconds, converged_count, iter_total)
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
    end subroutine print_asset_time_summary

    subroutine print_combo_summary(aic_wins, bic_wins, aic_rank_sum, bic_rank_sum, rank_count, fit_seconds)
        integer, intent(in) :: aic_wins(n_model,n_dist), bic_wins(n_model,n_dist)
        integer, intent(in) :: aic_rank_sum(n_model,n_dist), bic_rank_sum(n_model,n_dist)
        integer, intent(in) :: rank_count(n_model,n_dist)
        real(dp), intent(in) :: fit_seconds(n_model,n_dist)
        real(dp) :: mean_aic_rank, mean_bic_rank
        real(dp) :: freq
        integer :: im, id, model_total, dist_total, grand_total

        write(*,'(/,A)') "Model-distribution selection summary"
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
                if (model_total > 0) then
                    freq = real(aic_wins(im,id), dp) / real(model_total, dp)
                else
                    freq = 0.0_dp
                end if
                write(*,'(A16,1X,A10,1X,I9,1X,F9.3)') &
                    trim(models(im)), trim(dists(id)), aic_wins(im,id), freq
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
            if (grand_total > 0) then
                freq = real(dist_total, dp) / real(grand_total, dp)
            else
                freq = 0.0_dp
            end if
            write(*,'(A10,1X,I9,1X,F9.3)') trim(dists(id)), dist_total, freq
        end do
        write(*,'(A)') repeat("-", 42)
    end subroutine print_combo_summary

end program xfit_gen_garch_dist_returns

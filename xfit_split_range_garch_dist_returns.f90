! Fit split close-to-open/open-to-close GARCH models with multiple innovation
! distributions and high-low range information to adjusted OHLC prices.
!
! Usage:
!   xfit_split_range_garch_dist_returns.exe [ohlc_prices_file] [output_csv]
!
! For each asset, CO_t = log(Open_t / Close_{t-1}), OC_t = log(Close_t / Open_t),
! and CC_t = log(Close_t / Close_{t-1}) are demeaned.  The daily GARCH variance
! path is driven by CC_t, and co_frac allocates the daily variance between CO
! and OC.  CO and OC standardized innovations are treated as conditionally
! independent draws from the same selected innovation distribution.  The
! high-low range contributes a Parkinson-style quasi-likelihood term tied to
! the OC variance.

module split_range_garch_dist_returns_mod
    use kind_mod, only: dp
    use csv_mod, only: read_ohlc_csv, print_price_sample_info
    use stats_mod, only: mean
    use garch_fit_dist_mod, only: dist_param_count, dist_nu_value, dist_xi_value, dist_alpha_value
    use garch_split_fit_dist_mod, only: split_garch_dist_fit_result_t, fit_split_garch_range_dist_model, &
        split_model_param_count, split_garch_vol_forecast
    implicit none

    character(len=*), parameter :: default_prices_file = "prices_ohlc.csv"
    integer, parameter :: max_iter = 120
    real(dp), parameter :: gtol = 1.0e-6_dp
    character(len=16), parameter :: models(*) = [character(len=16) :: &
        "SPLIT_SYMM", "SPLIT_NAGARCH", "SPLIT_GJR", "SPLIT_FGTWIST"]
    character(len=10), parameter :: dists(*) = [character(len=10) :: &
        "NORMAL", "T", "SECH", "GED", "LAPLACE", "LOGISTIC", "NIG", "FS_SKEWT"]
    integer, parameter :: n_model = size(models)
    integer, parameter :: n_dist = size(dists)
    integer, parameter :: n_combo = n_model*n_dist

contains

    subroutine run_split_range_garch_dist_returns()
    integer, allocatable :: dates(:)
    character(len=32), allocatable :: asset_names(:)
    real(dp), allocatable :: open_prices(:,:), close_prices(:,:), high_prices(:,:), low_prices(:,:)
    real(dp), allocatable :: ret_cc(:), ret_co(:), ret_oc(:), range_var(:), vol(:)
    type(split_garch_dist_fit_result_t) :: rows(n_combo)
    integer :: row_model_idx(n_combo), row_dist_idx(n_combo), row_aic_rank(n_combo), row_bic_rank(n_combo)
    integer :: aic_wins(n_model,n_dist), bic_wins(n_model,n_dist)
    integer :: nprices, nassets, nobs, iasset, imod, idist, ifit, jfit, nfit
    integer :: nparam, clock_start, clock_end, clock_rate, fit_clock_start, fit_clock_end, csv_unit
    real(dp) :: row_logl(n_combo), row_aic(n_combo), row_bic(n_combo), row_fit_sec(n_combo), row_vol_ann(n_combo)
    real(dp) :: elapsed_s
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

    call read_ohlc_csv(trim(prices_file), dates, asset_names, open_prices, close_prices, &
                       high_prices=high_prices, low_prices=low_prices)
    nprices = size(close_prices, 1)
    nassets = size(close_prices, 2)
    nobs = nprices - 1
    allocate(ret_cc(nobs), ret_co(nobs), ret_oc(nobs), range_var(nobs), vol(nobs))

    call print_price_sample_info(trim(prices_file), dates, nassets)
    write(*,'(A)') "Split CO/OC/range likelihood: CO and OC distribution terms plus a Parkinson high-low range quasi-likelihood."
    write(*,'(A)') "h_co = co_frac*h_daily, h_oc = (1-co_frac)*h_daily; news impact uses CC returns."
    if (write_csv) then
        open(newunit=csv_unit, file=trim(output_csv), status="replace", action="write")
        write(csv_unit,'(A)') "model,dist,asset,omega,alpha,gamma,beta,theta,twist,co_frac,nu,xi,dist_alpha,persist,vol_ann_pct,logL,AIC,BIC,iter,conv,AIC_rank,BIC_rank,fit_sec"
    end if

    write(*,'(A)') "Model            Dist       Asset        omega   alpha   gamma    beta   theta   twist co_frac      nu      xi  dist_alpha  persist  vol_ann%        logL         AIC         BIC  iter conv AIC_rank BIC_rank    fit_sec"
    write(*,'(A)') repeat("-", 214)

    do iasset = 1, nassets
        ret_cc = log(close_prices(2:nprices,iasset) / close_prices(1:nprices-1,iasset))
        ret_cc = ret_cc - mean(ret_cc)
        ret_co = log(open_prices(2:nprices,iasset) / close_prices(1:nprices-1,iasset))
        ret_co = ret_co - mean(ret_co)
        ret_oc = log(close_prices(2:nprices,iasset) / open_prices(2:nprices,iasset))
        ret_oc = ret_oc - mean(ret_oc)
        range_var = (log(high_prices(2:nprices,iasset) / low_prices(2:nprices,iasset)))**2 / (4.0_dp*log(2.0_dp))
        nfit = 0

        do imod = 1, n_model
            do idist = 1, n_dist
                nfit = nfit + 1
                row_model_idx(nfit) = imod
                row_dist_idx(nfit) = idist
                call system_clock(fit_clock_start)
                call fit_split_garch_range_dist_model(trim(models(imod)), trim(dists(idist)), ret_cc, ret_co, &
                                                      ret_oc, range_var, max_iter, gtol, rows(nfit))
                call system_clock(fit_clock_end)
                row_fit_sec(nfit) = real(fit_clock_end - fit_clock_start, dp) / real(clock_rate, dp)
                row_logl(nfit) = -real(nobs, dp)*rows(nfit)%nll
                nparam = split_model_param_count(trim(models(imod))) + dist_param_count(trim(dists(idist)))
                row_aic(nfit) = 2.0_dp*real(nparam, dp) - 2.0_dp*row_logl(nfit)
                row_bic(nfit) = log(real(nobs, dp))*real(nparam, dp) - 2.0_dp*row_logl(nfit)
                call split_garch_vol_forecast(ret_cc, trim(models(imod)), rows(nfit)%params, rows(nfit)%persist, vol)
                row_vol_ann(nfit) = sqrt(sum(vol**2) / real(nobs, dp))
            end do
        end do

        do ifit = 1, nfit
            row_aic_rank(ifit) = 1
            row_bic_rank(ifit) = 1
            do jfit = 1, nfit
                if (row_aic(jfit) < row_aic(ifit)) row_aic_rank(ifit) = row_aic_rank(ifit) + 1
                if (row_bic(jfit) < row_bic(ifit)) row_bic_rank(ifit) = row_bic_rank(ifit) + 1
            end do
            if (row_aic_rank(ifit) == 1) aic_wins(row_model_idx(ifit), row_dist_idx(ifit)) = &
                aic_wins(row_model_idx(ifit), row_dist_idx(ifit)) + 1
            if (row_bic_rank(ifit) == 1) bic_wins(row_model_idx(ifit), row_dist_idx(ifit)) = &
                bic_wins(row_model_idx(ifit), row_dist_idx(ifit)) + 1
        end do

        do ifit = 1, nfit
            call print_fit_row(rows(ifit), trim(asset_names(iasset)), row_vol_ann(ifit), row_logl(ifit), &
                               row_aic(ifit), row_bic(ifit), row_aic_rank(ifit), row_bic_rank(ifit), &
                               row_fit_sec(ifit), write_csv, csv_unit)
        end do
    end do

    if (write_csv) then
        close(csv_unit)
        write(*,'(/,A,A)') "Wrote fit table CSV: ", trim(output_csv)
    end if
    call print_selection_summary(aic_wins, bic_wins)
    call system_clock(clock_end)
    elapsed_s = real(clock_end - clock_start, dp) / real(clock_rate, dp)
    write(*,'(/,A,F10.3,A)') "Elapsed wall time: ", elapsed_s, " seconds"
    end subroutine run_split_range_garch_dist_returns

    subroutine print_fit_row(row, asset, vol_ann, logl, aic, bic, aic_rank, bic_rank, fit_sec, write_csv, csv_unit)
        type(split_garch_dist_fit_result_t), intent(in) :: row
        character(len=*), intent(in) :: asset
        real(dp), intent(in) :: vol_ann, logl, aic, bic, fit_sec
        integer, intent(in) :: aic_rank, bic_rank, csv_unit
        logical, intent(in) :: write_csv

        write(*,'(A16,1X,A10,1X,A9,ES12.3,6F8.4,1X,A8,1X,A8,1X,A10,F9.4,F10.2,F12.2,F12.2,F12.2,I6,1X,L1,2I9,1X,F10.3)') &
            trim(row%model), trim(row%dist), trim(asset), row%params%omega, row%params%alpha, row%params%gamma, &
            row%params%beta, row%params%theta, row%params%twist, row%co_frac, &
            dist_nu_value(trim(row%dist), row%shape), dist_xi_value(trim(row%dist), row%shape2), &
            dist_alpha_value(trim(row%dist), row%shape), row%persist, vol_ann, logl, aic, bic, row%niter, &
            row%converged, aic_rank, bic_rank, fit_sec
        if (write_csv) then
            write(csv_unit,'(A,",",A,",",A,",",ES16.8,",",F12.6,",",F12.6,",",F12.6,",",F12.6,",",F12.6,",",F12.6,",",A,",",A,",",A,",",F12.6,",",F12.6,",",F16.6,",",F16.6,",",F16.6,",",I0,",",L1,",",I0,",",I0,",",F12.6)') &
                trim(row%model), trim(row%dist), trim(asset), row%params%omega, row%params%alpha, row%params%gamma, &
                row%params%beta, row%params%theta, row%params%twist, row%co_frac, &
                dist_nu_value(trim(row%dist), row%shape), dist_xi_value(trim(row%dist), row%shape2), &
                dist_alpha_value(trim(row%dist), row%shape), row%persist, vol_ann, logl, aic, bic, row%niter, &
                row%converged, aic_rank, bic_rank, fit_sec
        end if
    end subroutine print_fit_row

    subroutine print_selection_summary(aic_wins, bic_wins)
        integer, intent(in) :: aic_wins(n_model,n_dist), bic_wins(n_model,n_dist)
        integer :: im, id

        write(*,'(/,A)') "Split CO/OC/range AIC/BIC win counts"
        write(*,'(A)') repeat("-", 48)
        write(*,'(A16,1X,A10,1X,A8,1X,A8)') "Model", "Dist", "AIC_wins", "BIC_wins"
        write(*,'(A)') repeat("-", 48)
        do im = 1, n_model
            do id = 1, n_dist
                write(*,'(A16,1X,A10,1X,I8,1X,I8)') trim(models(im)), trim(dists(id)), aic_wins(im,id), bic_wins(im,id)
            end do
        end do
        write(*,'(A)') repeat("-", 48)
    end subroutine print_selection_summary

end module split_range_garch_dist_returns_mod

program xfit_split_range_garch_dist_returns
    use split_range_garch_dist_returns_mod, only: run_split_range_garch_dist_returns
    implicit none

    call run_split_range_garch_dist_returns()
end program xfit_split_range_garch_dist_returns

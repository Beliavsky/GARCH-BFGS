! Fit the return-volatility models that are directly comparable to rugarch
! built-in GARCH-family models.

program xfit_rugarch_builtin_symm_returns
    use kind_mod,        only: dp
    use csv_mod,         only: read_price_csv, print_price_sample_info
    use stats_mod,       only: mean
    use garch_types_mod, only: garch_params_t, garch_fit_result_t
    use garch_forecast_mod, only: finalize_return_garch_fit
    use model_selection_mod, only: print_model_selection_counts, print_model_fit_times
    use garch_fit_mod,   only: fit_symm_garch, fit_symm_garch_pq, fit_figarch, fit_csgarch, fit_tgarch, fit_avgarch
    implicit none

    character(len=*), parameter :: prices_file = "spy_efa_eem_tlt_lqd.csv"
    character(len=*), parameter :: csv_file = "rugarch/fortran_rugarch_builtin_symm_results.csv"
    real(dp), parameter :: trading_days = 252.0_dp
    integer,  parameter :: max_iter = 500
    real(dp), parameter :: gtol = 1.0e-7_dp
    character(len=16), parameter :: models(*) = [character(len=16) :: &
        "SYMM_GARCH", "SYMM_GARCH_2_1", "SYMM_GARCH_1_2", "SYMM_GARCH_2_2", &
        "FIGARCH", "CSGARCH", "TGARCH", "AVGARCH"]
    integer, parameter :: n_model = size(models)

    type row_t
        character(len=16) :: model = ""
        character(len=32) :: asset = ""
        type(garch_fit_result_t) :: fit
        real(dp) :: seconds = 0.0_dp
        integer :: aic_rank = 0
        integer :: bic_rank = 0
    end type row_t

    integer, allocatable :: dates(:)
    character(len=32), allocatable :: col_names(:)
    real(dp), allocatable :: prices(:,:), ret(:)
    type(row_t), allocatable :: rows(:)
    integer :: nprices, ncols, nobs, icol, imod, irow
    integer :: clock_start, clock_end, clock_rate
    real(dp) :: elapsed_s, ret_mean
    integer :: aic_wins(n_model), bic_wins(n_model), param_count(n_model)
    integer :: aic_rank_sum(n_model), bic_rank_sum(n_model), rank_count(n_model)
    real(dp) :: fit_seconds(n_model)
    character(len=256) :: aic_symbols(n_model)

    call system_clock(clock_start, clock_rate)
    aic_wins = 0
    bic_wins = 0
    aic_rank_sum = 0
    bic_rank_sum = 0
    rank_count = 0
    param_count = 0
    fit_seconds = 0.0_dp
    aic_symbols = ""

    call read_price_csv(prices_file, dates, col_names, prices)
    nprices = size(prices, 1)
    ncols = size(prices, 2)
    nobs = nprices - 1
    allocate(ret(nobs), rows(ncols*n_model))

    irow = 0
    do icol = 1, ncols
        ret = log(prices(2:nprices,icol) / prices(1:nprices-1,icol))
        ret_mean = mean(ret)
        ret = ret - ret_mean
        do imod = 1, n_model
            irow = irow + 1
            call fit_model(models(imod), trim(col_names(icol)), ret, rows(irow))
            param_count(imod) = rows(irow)%fit%nparam
            fit_seconds(imod) = fit_seconds(imod) + rows(irow)%seconds
        end do
    end do

    call rank_results(rows, irow)
    call update_summary(rows, irow, aic_wins, bic_wins, aic_rank_sum, bic_rank_sum, rank_count, aic_symbols)

    call print_price_sample_info(prices_file, dates, ncols)
    write(*,'(A)') "Models: rugarch built-in GARCH family"
    call print_rows(rows, irow)
    call print_model_selection_counts(models, param_count, aic_wins, bic_wins, aic_symbols, &
                                      aic_rank_sum, bic_rank_sum, rank_count)
    call print_model_fit_times(models, param_count, fit_seconds)
    call write_results_csv(rows, irow)

    call system_clock(clock_end)
    elapsed_s = real(clock_end - clock_start, dp) / real(clock_rate, dp)
    write(*,'(/,A,A)') "Wrote CSV results: ", csv_file
    write(*,'(A,F10.3,A)') "Elapsed wall time: ", elapsed_s, " seconds"

contains

    subroutine fit_model(model_name, asset, y, row)
        character(len=*), intent(in) :: model_name, asset
        real(dp), intent(in) :: y(:)
        type(row_t), intent(out) :: row
        type(garch_params_t) :: params
        integer :: niter, t0, t1, rate
        real(dp) :: fopt
        logical :: converged

        call system_clock(t0, rate)
        row = row_t()
        row%model = model_name
        row%asset = asset

        select case (trim(model_name))
        case ("SYMM_GARCH")
            call fit_symm_garch(y, max_iter, gtol, fopt, params, niter, converged)
        case ("SYMM_GARCH_2_1")
            call fit_symm_garch_pq(y, 2, 1, max_iter, gtol, fopt, params, niter, converged)
        case ("SYMM_GARCH_1_2")
            call fit_symm_garch_pq(y, 1, 2, max_iter, gtol, fopt, params, niter, converged)
        case ("SYMM_GARCH_2_2")
            call fit_symm_garch_pq(y, 2, 2, max_iter, gtol, fopt, params, niter, converged)
        case ("FIGARCH")
            call fit_figarch(y, max_iter, gtol, fopt, params, niter, converged)
        case ("CSGARCH")
            call fit_csgarch(y, max_iter, gtol, fopt, params, niter, converged)
        case ("TGARCH")
            call fit_tgarch(y, max_iter, gtol, fopt, params, niter, converged)
        case ("AVGARCH")
            call fit_avgarch(y, max_iter, gtol, fopt, params, niter, converged)
        case default
            write(*,'(A,A)') "Unknown model: ", trim(model_name)
            error stop
        end select

        call finalize_return_garch_fit(model_name, y, params, fopt, niter, converged, trading_days, row%fit)
        call system_clock(t1)
        row%seconds = real(t1 - t0, dp) / real(rate, dp)
    end subroutine fit_model

    subroutine rank_results(rows, nrow)
        type(row_t), intent(inout) :: rows(:)
        integer, intent(in) :: nrow
        integer :: i, j

        do i = 1, nrow
            rows(i)%aic_rank = 1
            rows(i)%bic_rank = 1
            do j = 1, nrow
                if (trim(rows(j)%asset) /= trim(rows(i)%asset)) cycle
                if (rows(j)%fit%aic < rows(i)%fit%aic) rows(i)%aic_rank = rows(i)%aic_rank + 1
                if (rows(j)%fit%bic < rows(i)%fit%bic) rows(i)%bic_rank = rows(i)%bic_rank + 1
            end do
        end do
    end subroutine rank_results

    subroutine update_summary(rows, nrow, aic_wins, bic_wins, aic_rank_sum, bic_rank_sum, rank_count, aic_symbols)
        type(row_t), intent(in) :: rows(:)
        integer, intent(in) :: nrow
        integer, intent(inout) :: aic_wins(:), bic_wins(:), aic_rank_sum(:), bic_rank_sum(:), rank_count(:)
        character(len=*), intent(inout) :: aic_symbols(:)
        integer :: i, imod

        do i = 1, nrow
            imod = model_index(rows(i)%model)
            if (imod <= 0) cycle
            if (rows(i)%aic_rank == 1) then
                aic_wins(imod) = aic_wins(imod) + 1
                if (len_trim(aic_symbols(imod)) == 0) then
                    aic_symbols(imod) = trim(rows(i)%asset)
                else
                    aic_symbols(imod) = trim(aic_symbols(imod)) // " " // trim(rows(i)%asset)
                end if
            end if
            if (rows(i)%bic_rank == 1) bic_wins(imod) = bic_wins(imod) + 1
            aic_rank_sum(imod) = aic_rank_sum(imod) + rows(i)%aic_rank
            bic_rank_sum(imod) = bic_rank_sum(imod) + rows(i)%bic_rank
            rank_count(imod) = rank_count(imod) + 1
        end do
    end subroutine update_summary

    integer function model_index(model_name)
        character(len=*), intent(in) :: model_name
        integer :: i

        model_index = 0
        do i = 1, n_model
            if (trim(model_name) == trim(models(i))) then
                model_index = i
                return
            end if
        end do
    end function model_index

    subroutine print_rows(rows, nrow)
        type(row_t), intent(in) :: rows(:)
        integer, intent(in) :: nrow
        integer :: i

        write(*,'(A)') "Model            Asset        omega   alpha   gamma    beta   theta   twist  persist  vol_ann%        logL         AIC         BIC  iter conv    skew   ekurt AIC_rank BIC_rank"
        write(*,'(A)') repeat("-", 172)
        do i = 1, nrow
            write(*,'(A16,1X,A9,ES12.3,6F8.4,F10.2,F12.2,F12.2,F12.2,I6,1X,L1,2F8.3,2I9)') &
                trim(rows(i)%model), trim(rows(i)%asset), rows(i)%fit%params%omega, rows(i)%fit%params%alpha, &
                rows(i)%fit%params%gamma, rows(i)%fit%params%beta, rows(i)%fit%params%theta, &
                rows(i)%fit%params%twist, rows(i)%fit%persist, rows(i)%fit%vol_ann, rows(i)%fit%logl, &
                rows(i)%fit%aic, rows(i)%fit%bic, rows(i)%fit%niter, rows(i)%fit%converged, &
                rows(i)%fit%skew, rows(i)%fit%ekurt, rows(i)%aic_rank, rows(i)%bic_rank
        end do
    end subroutine print_rows

    subroutine write_results_csv(rows, nrow)
        type(row_t), intent(in) :: rows(:)
        integer, intent(in) :: nrow
        integer :: unit, i

        open(newunit=unit, file=csv_file, status="replace", action="write")
        write(unit,'(A)') "model,omega,alpha,gamma,beta,theta,twist,persist,vol_ann,logl,aic,bic,nparam,niter,conv,skew,ekurt,sec,asset,aic_rank,bic_rank"
        do i = 1, nrow
            write(unit,'(A,",",ES24.16,",",ES24.16,",",ES24.16,",",ES24.16,",",ES24.16,",",ES24.16,",",ES24.16,",",ES24.16,",",ES24.16,",",ES24.16,",",ES24.16,",",I0,",",I0,",",L1,",",ES24.16,",",ES24.16,",",ES24.16,",",A,",",I0,",",I0)') &
                trim(rows(i)%model), rows(i)%fit%params%omega, rows(i)%fit%params%alpha, rows(i)%fit%params%gamma, &
                rows(i)%fit%params%beta, rows(i)%fit%params%theta, rows(i)%fit%params%twist, rows(i)%fit%persist, &
                rows(i)%fit%vol_ann, rows(i)%fit%logl, rows(i)%fit%aic, rows(i)%fit%bic, rows(i)%fit%nparam, &
                rows(i)%fit%niter, rows(i)%fit%converged, rows(i)%fit%skew, rows(i)%fit%ekurt, rows(i)%seconds, &
                trim(rows(i)%asset), rows(i)%aic_rank, rows(i)%bic_rank
        end do
        close(unit)
    end subroutine write_results_csv

end program xfit_rugarch_builtin_symm_returns

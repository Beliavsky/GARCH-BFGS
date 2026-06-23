! Fit GARCH-family models × innovation distributions using a rolling training window.
! Each window is warm-started from the previous window's parameter estimates.
! Models are ranked by cumulative out-of-sample log-likelihood.
!
! Rolling scheme:
!   - Fit on ret(t0 : t0+nroll-1) (length nroll).
!   - Evaluate OOS log-likelihood on ret(t0+nroll : t0+nroll+nstep-1) (length nstep).
!   - Advance t0 by nstep and repeat.
!
! Usage: xfit_gen_garch_roll.exe [prices_file]

program xfit_gen_garch_roll
    use date_mod,           only: print_program_header, date_label
    use kind_mod,           only: dp
    use csv_mod,            only: read_price_csv, print_price_sample_info
    use stats_mod,          only: mean
    use garch_fit_dist_mod, only: garch_dist_fit_result_t, fit_garch_dist_model, &
                                  garch_dist_oos_nll
    implicit none

    character(len=*), parameter :: default_prices_file = "spy_efa_eem_tlt_lqd.csv"
    integer,  parameter :: max_assets = 3, max_iter = 100
    real(dp), parameter :: gtol = 1.0e-7_dp

    ! ---- user settings ----
    integer, parameter :: nroll = 1000   ! training window length (observations)
    integer, parameter :: nstep = 252     ! periods between refits; also OOS window length
    character(len=16), parameter :: models(*) = [character(len=16) :: &
        "SYMM_GARCH", "NAGARCH"]
    character(len=10), parameter :: dists(*) = [character(len=10) :: &
        "NORMAL", "T"] ! , "FS_SKEWT"]
    ! -----------------------

    integer, parameter :: n_model = size(models)
    integer, parameter :: n_dist  = size(dists)
    integer, parameter :: n_combo = n_model * n_dist

    integer,          allocatable :: dates(:)
    character(len=32),allocatable :: col_names(:)
    real(dp),         allocatable :: prices(:,:), ret(:)
    type(garch_dist_fit_result_t), allocatable :: cur_res(:), prev_res(:)
    real(dp), allocatable :: oos_logl_total(:)
    logical,  allocatable :: prev_fitted(:)
    real(dp) :: oos_logl(n_combo), oos_nll_val
    integer  :: rank_order(n_combo)
    integer  :: row_model_idx(n_combo), row_dist_idx(n_combo)
    integer  :: nprices, ncols, nobs
    integer  :: icol, imod, idist, ifit, i, j, tmp_idx
    integer  :: t0, t1, t2, nwindows, nobs_oos_total
    integer  :: clock_start, clock_end, clock_rate
    real(dp) :: ret_mean, elapsed_s
    character(len=256) :: prices_file

    call print_program_header("xfit_gen_garch_roll.f90")
    prices_file = default_prices_file
    if (command_argument_count() >= 1) call get_command_argument(1, prices_file)

    call system_clock(clock_start, clock_rate)
    call read_price_csv(prices_file, dates, col_names, prices, max_col=max_assets)
    nprices = size(prices, 1)
    ncols   = size(prices, 2)
    nobs    = nprices - 1
    allocate(ret(nobs))
    allocate(cur_res(n_combo), prev_res(n_combo))
    allocate(oos_logl_total(n_combo), prev_fitted(n_combo))

    ifit = 0
    do imod = 1, n_model
        do idist = 1, n_dist
            ifit = ifit + 1
            row_model_idx(ifit) = imod
            row_dist_idx(ifit)  = idist
        end do
    end do

    call print_price_sample_info(prices_file, dates, ncols)
    write(*,'(/,A,I0,A,I0,A)') &
        "Rolling: nroll=", nroll, "  nstep=", nstep, "  (warm starts)"

    do icol = 1, ncols
        ret      = log(prices(2:nprices,icol) / prices(1:nprices-1,icol))
        ret_mean = mean(ret)
        ret      = ret - ret_mean

        oos_logl_total = 0.0_dp
        prev_fitted    = .false.
        nwindows       = 0
        nobs_oos_total = 0

        if (nobs < nroll + nstep) then
            write(*,'(/,A,A,A)') "Skipping ", trim(col_names(icol)), &
                ": nobs < nroll + nstep"
            cycle
        end if

        write(*,'(/,A,A)') "=== Asset: ", trim(col_names(icol))
        call print_window_header()

        t0 = 1
        do while (t0 + nroll + nstep - 1 <= nobs)
            t1 = t0 + nroll - 1
            t2 = t0 + nroll + nstep - 1
            nwindows       = nwindows + 1
            nobs_oos_total = nobs_oos_total + nstep

            do ifit = 1, n_combo
                imod  = row_model_idx(ifit)
                idist = row_dist_idx(ifit)
                if (prev_fitted(ifit)) then
                    call fit_garch_dist_model(models(imod), dists(idist), &
                        ret(t0:t1), max_iter, gtol, cur_res(ifit), &
                        start_params = prev_res(ifit)%params, &
                        start_shape  = prev_res(ifit)%shape, &
                        start_shape2 = prev_res(ifit)%shape2)
                else
                    call fit_garch_dist_model(models(imod), dists(idist), &
                        ret(t0:t1), max_iter, gtol, cur_res(ifit))
                end if
                oos_nll_val = garch_dist_oos_nll(models(imod), dists(idist), &
                    ret(t0:t2), nroll, nstep, cur_res(ifit))
                oos_logl(ifit)        = -real(nstep, dp) * oos_nll_val
                oos_logl_total(ifit)  = oos_logl_total(ifit) + oos_logl(ifit)
                prev_res(ifit)        = cur_res(ifit)
                prev_fitted(ifit)     = cur_res(ifit)%converged
            end do

            call print_window_row(nwindows, t0, t1, t2, oos_logl)
            t0 = t0 + nstep
        end do

        ! sort indices by descending total OOS logL
        do ifit = 1, n_combo
            rank_order(ifit) = ifit
        end do
        do i = 1, n_combo - 1
            do j = i + 1, n_combo
                if (oos_logl_total(rank_order(j)) > oos_logl_total(rank_order(i))) then
                    tmp_idx        = rank_order(i)
                    rank_order(i)  = rank_order(j)
                    rank_order(j)  = tmp_idx
                end if
            end do
        end do

        call print_oos_summary(col_names(icol), nwindows, nobs_oos_total, &
            oos_logl_total, rank_order)
    end do

    call system_clock(clock_end)
    elapsed_s = real(clock_end - clock_start, dp) / real(clock_rate, dp)
    write(*,'(/,A,F10.3,A)') "Elapsed wall time:", elapsed_s, " seconds"

    deallocate(dates, col_names, prices, ret)
    deallocate(cur_res, prev_res, oos_logl_total, prev_fitted)

contains

    pure function model_abbrev(im) result(s)
        integer, intent(in) :: im
        character(len=4) :: s
        select case (trim(models(im)))
        case ("SYMM_GARCH", "GARCH"); s = "SG"
        case ("NAGARCH");              s = "NAG"
        case ("GJR_GARCH", "GJR");    s = "GJR"
        case ("GJR_SIGNED");           s = "GJRS"
        case ("EGARCH");               s = "EG"
        case ("QGARCH");               s = "QG"
        case default;                  s = models(im)(1:4)
        end select
    end function model_abbrev

    pure function dist_abbrev(id) result(s)
        integer, intent(in) :: id
        character(len=6) :: s
        select case (trim(dists(id)))
        case ("NORMAL");   s = "N"
        case ("T");        s = "T"
        case ("FS_SKEWT"); s = "FS"
        case ("SECH");     s = "SECH"
        case ("GED");      s = "GED"
        case ("LAPLACE");  s = "LAP"
        case ("LOGISTIC"); s = "LOG"
        case ("NIG");      s = "NIG"
        case default;      s = dists(id)(1:min(6, len_trim(dists(id))))
        end select
    end function dist_abbrev

    pure function combo_label(im, id) result(lbl)
        integer, intent(in) :: im, id
        character(len=10) :: lbl
        lbl = trim(model_abbrev(im)) // "/" // trim(dist_abbrev(id))
    end function combo_label

    subroutine print_window_header()
        integer :: ic, w
        ! width: 5 (Win) + 4*11 (dates) + n_combo*10 (OOS cols)
        w = 5 + 4*11 + n_combo*10
        write(*,'(A)') repeat("-", w)
        write(*,'(A5,4(1X,A10))', advance='no') &
            "Win", "TrnStart", "TrnEnd", "OOS_Start", "OOS_End"
        do ic = 1, n_combo
            write(*,'(1X,A9)', advance='no') &
                combo_label(row_model_idx(ic), row_dist_idx(ic))
        end do
        write(*,*)
        write(*,'(A)') repeat("-", w)
    end subroutine print_window_header

    subroutine print_window_row(iwin, t0, t1, t2, oos_logl)
        integer,  intent(in) :: iwin, t0, t1, t2
        real(dp), intent(in) :: oos_logl(n_combo)
        integer :: ic
        ! dates(t+1) is the closing-price date for ret(t)
        write(*,'(I5,1X,A10,1X,A10,1X,A10,1X,A10)', advance='no') &
            iwin, &
            date_label(dates(t0 + 1)), &
            date_label(dates(t1 + 1)), &
            date_label(dates(t1 + 2)), &
            date_label(dates(t2 + 1))
        do ic = 1, n_combo
            write(*,'(1X,F9.4)', advance='no') oos_logl(ic) / real(nstep, dp)
        end do
        write(*,*)
    end subroutine print_window_row

    subroutine print_oos_summary(asset, nwin, nobs_oos, total_logl, rank_order)
        character(len=*), intent(in) :: asset
        integer,  intent(in) :: nwin, nobs_oos, rank_order(n_combo)
        real(dp), intent(in) :: total_logl(n_combo)
        integer :: k, ifit, im, id
        write(*,'(/,A,1X,A,A,I0,A,I0,A)') &
            "OOS log-likelihood summary:", trim(asset), &
            "  (", nwin, " windows,", nobs_oos, " OOS obs)"
        write(*,'(A)') repeat("-", 66)
        write(*,'(A16,1X,A10,1X,A14,1X,A10,1X,A4)') &
            "Model", "Dist", "Total_OOS_logL", "Per_obs", "Rank"
        write(*,'(A)') repeat("-", 66)
        do k = 1, n_combo
            ifit = rank_order(k)
            im   = row_model_idx(ifit)
            id   = row_dist_idx(ifit)
            write(*,'(A16,1X,A10,1X,F14.2,1X,F10.5,1X,I4)') &
                trim(models(im)), trim(dists(id)), &
                total_logl(ifit), &
                total_logl(ifit) / real(max(nobs_oos, 1), dp), &
                k
        end do
        write(*,'(A)') repeat("-", 66)
    end subroutine print_oos_summary

end program xfit_gen_garch_roll

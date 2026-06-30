! Find the rolling-window length (nroll = k*nstep) that maximises per-observation
! out-of-sample log-likelihood for GARCH-family × innovation-distribution combinations.
!
! For each candidate nroll = k*nstep (k = nroll_min_mult..nroll_max_mult):
!   - Fit on ret(t0 : t0+nroll-1), evaluate OOS on ret(t0+nroll : t0+nroll+nstep-1)
!   - Advance t0 by nstep; warm-start from the previous window's parameter estimates
!   - Accumulate per-observation OOS log-likelihood
!
! Summary table: rows = nroll values, columns = model/dist combos.
! Best nroll (highest per-obs OOS logL) identified per combo.
!
! Usage: xfit_garch_roll_period.exe [prices_file]

program xfit_garch_roll_period
    use date_mod,           only: print_program_header
    use kind_mod,           only: dp
    use csv_mod,            only: read_price_csv, print_price_sample_info
    use stats_mod,          only: mean
    use garch_fit_dist_mod, only: garch_dist_fit_result_t, fit_garch_dist_model, &
                                  garch_dist_oos_nll, model_abbrev, dist_abbrev
    implicit none

    character(len=*), parameter :: default_prices_file = "spy_efa_eem_tlt_lqd.csv"
    integer,  parameter :: max_assets = 10000, max_iter = 100
    real(dp), parameter :: gtol = 1.0e-7_dp

    ! ---- user settings ----
    integer, parameter :: nstep          = 252   ! OOS evaluation window; also step size
    integer, parameter :: nroll_min_mult =  3    ! smallest nroll = nroll_min_mult * nstep
    integer, parameter :: nroll_max_mult =  10   ! largest  nroll = nroll_max_mult * nstep
    character(len=16), parameter :: models(*) = [character(len=16) :: &
        "SYMM_GARCH", "NAGARCH"]
    character(len=10), parameter :: dists(*) = [character(len=10) :: &
        "NORMAL", "T"]
    ! -----------------------

    integer, parameter :: n_model = size(models)
    integer, parameter :: n_dist  = size(dists)
    integer, parameter :: n_combo = n_model * n_dist
    integer, parameter :: n_mult  = nroll_max_mult - nroll_min_mult + 1

    integer,           allocatable :: dates(:)
    character(len=32), allocatable :: col_names(:)
    real(dp),          allocatable :: prices(:,:), ret(:)
    type(garch_dist_fit_result_t), allocatable :: cur_res(:), prev_res(:)
    real(dp), allocatable :: oos_logl_total(:)
    logical,  allocatable :: prev_fitted(:)

    real(dp) :: oos_per_obs(n_mult, n_combo)
    integer  :: nobs_oos_by_mult(n_mult), nwin_by_mult(n_mult)
    integer  :: row_model_idx(n_combo), row_dist_idx(n_combo)

    integer  :: nprices, ncols, nobs
    integer  :: icol, imod, idist, ifit, imult, mult, nroll
    integer  :: t0, t1, t2, nwindows, nobs_oos_total
    integer  :: clock_start, clock_end, clock_rate
    real(dp) :: ret_mean, elapsed_s, oos_nll_val
    character(len=256) :: prices_file

    call print_program_header("xfit_garch_roll_period.f90")
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
    write(*,'(A,I0,A,I0,A,I0)') &
        "nstep=", nstep, "  nroll range: ", nroll_min_mult*nstep, &
        " to ", nroll_max_mult*nstep

    do icol = 1, ncols
        ret      = log(prices(2:nprices,icol) / prices(1:nprices-1,icol))
        ret_mean = mean(ret)
        ret      = ret - ret_mean

        write(*,'(/,A,A)') "=== Asset: ", trim(col_names(icol))

        do imult = 1, n_mult
            mult  = nroll_min_mult + imult - 1
            nroll = mult * nstep

            if (nobs < nroll + nstep) then
                write(*,'(2X,A,I0,A)') "nroll=", nroll, &
                    ": skipped (insufficient observations)"
                nobs_oos_by_mult(imult) = 0
                nwin_by_mult(imult)     = 0
                oos_per_obs(imult, :)   = -huge(1.0_dp)
                cycle
            end if

            oos_logl_total = 0.0_dp
            prev_fitted    = .false.
            nwindows       = 0
            nobs_oos_total = 0

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
                    oos_logl_total(ifit) = oos_logl_total(ifit) &
                                         - real(nstep, dp) * oos_nll_val
                    prev_res(ifit)    = cur_res(ifit)
                    prev_fitted(ifit) = cur_res(ifit)%converged
                end do
                t0 = t0 + nstep
            end do

            nobs_oos_by_mult(imult) = nobs_oos_total
            nwin_by_mult(imult)     = nwindows
            do ifit = 1, n_combo
                oos_per_obs(imult, ifit) = &
                    oos_logl_total(ifit) / real(nobs_oos_total, dp)
            end do
            write(*,'(2X,A,I0,A,I0,A,I0,A)') &
                "nroll=", nroll, "  windows=", nwindows, &
                "  OOS obs=", nobs_oos_total, " done"
        end do

        call print_roll_summary(col_names(icol))
    end do

    call system_clock(clock_end)
    elapsed_s = real(clock_end - clock_start, dp) / real(clock_rate, dp)
    write(*,'(/,A,F10.3,A)') "Elapsed wall time:", elapsed_s, " seconds"

    deallocate(dates, col_names, prices, ret)
    deallocate(cur_res, prev_res, oos_logl_total, prev_fitted)

contains

    subroutine print_roll_summary(asset)
        character(len=*), intent(in) :: asset
        integer  :: imult_idx, ifit, im, id, mult_val, best_mult
        real(dp) :: best_val
        ! column layout: nroll=I7, nwin=I5, nOOS=I7 (19 chars), each combo=1X+F9.4 (10 chars)
        integer, parameter :: lhs_w = 19, combo_w = 10

        write(*,'(/,A,A)') "Per-obs OOS log-likelihood by roll period: ", trim(asset)
        write(*,'(A7,A5,A7)', advance='no') "nroll", "nwin", "nOOS"
        do ifit = 1, n_combo
            im = row_model_idx(ifit)
            id = row_dist_idx(ifit)
            write(*,'(1X,A9)', advance='no') combo_label(im, id)
        end do
        write(*,*)
        write(*,'(A)') repeat("-", lhs_w + n_combo * combo_w)

        do imult_idx = 1, n_mult
            mult_val = (nroll_min_mult + imult_idx - 1) * nstep
            write(*,'(I7,I5,I7)', advance='no') &
                mult_val, nwin_by_mult(imult_idx), nobs_oos_by_mult(imult_idx)
            do ifit = 1, n_combo
                if (nobs_oos_by_mult(imult_idx) > 0) then
                    write(*,'(1X,F9.4)', advance='no') oos_per_obs(imult_idx, ifit)
                else
                    write(*,'(1X,A9)', advance='no') "   --    "
                end if
            end do
            write(*,*)
        end do

        write(*,'(A)') repeat("-", lhs_w + n_combo * combo_w)
        write(*,'(A7,A5,A7)', advance='no') "Best", "", ""
        do ifit = 1, n_combo
            best_val  = -huge(1.0_dp)
            best_mult = 0
            do imult_idx = 1, n_mult
                if (nobs_oos_by_mult(imult_idx) > 0 .and. &
                    oos_per_obs(imult_idx, ifit) > best_val) then
                    best_val  = oos_per_obs(imult_idx, ifit)
                    best_mult = (nroll_min_mult + imult_idx - 1) * nstep
                end if
            end do
            write(*,'(1X,I9)', advance='no') best_mult
        end do
        write(*,*)
        write(*,'(A)') repeat("-", lhs_w + n_combo * combo_w)
    end subroutine print_roll_summary

    pure function combo_label(im, id) result(lbl)
        integer, intent(in) :: im, id
        character(len=10) :: lbl
        lbl = trim(model_abbrev(models(im))) // "/" // trim(dist_abbrev(dists(id)))
    end function combo_label

end program xfit_garch_roll_period

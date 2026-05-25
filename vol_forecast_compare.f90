! Utilities for comparing model volatility forecasts with implied volatility data.

module vol_forecast_compare_mod
    use kind_mod,    only: dp
    use strings_mod, only: uppercase
    use csv_mod,     only: date_label
    use time_series_compare_mod, only: transform_series, lag_series_one, aligned_corr, is_change_transform
    implicit none
    private
    public :: print_implied_vol_correlations

contains

    subroutine print_implied_vol_correlations(implied_vol_file, asset_names, model_names, fit_model_have, &
                                             forecast_dates, vol_forecasts, iv_dates, iv_names, iv_values, &
                                             map_assets, map_iv_indices, iv_transform, horizons, &
                                             extra_names, extra_dates, extra_values, extra_signs, &
                                             standardize_extra_by_prior_iv)
        character(len=*), intent(in) :: implied_vol_file
        character(len=*), intent(in) :: asset_names(:), model_names(:), iv_names(:)
        logical, intent(in) :: fit_model_have(:)
        integer, intent(in) :: forecast_dates(:), iv_dates(:)
        real(dp), intent(in) :: vol_forecasts(:,:,:), iv_values(:,:)
        character(len=*), intent(in) :: map_assets(:), map_iv_indices(:)
        character(len=*), optional, intent(in) :: iv_transform
        integer, optional, intent(in) :: horizons(:)
        character(len=*), optional, intent(in) :: extra_names(:)
        integer, optional, intent(in) :: extra_dates(:)
        real(dp), optional, intent(in) :: extra_values(:,:,:)
        real(dp), optional, intent(in) :: extra_signs(:)
        logical, optional, intent(in) :: standardize_extra_by_prior_iv

        integer :: imap, iasset, iiv, imodel, ih, iextra, nmatch, first_index, last_index
        integer :: max_rows, nrow, irow
        integer :: hvals(1)
        integer, allocatable :: use_horizons(:), iv_corr_dates(:), iv_lag_dates(:), garch_corr_dates(:), std_dates(:)
        integer, allocatable :: row_nmatch(:), row_first_index(:), row_last_index(:)
        real(dp), allocatable :: iv_corr_values(:), iv_lag_values(:), garch_corr_values(:), std_values(:)
        real(dp), allocatable :: row_corr(:)
        real(dp) :: corr, extra_sign
        logical, allocatable :: row_valid(:)
        logical :: do_standardize_extra
        character(len=16), allocatable :: row_model(:)
        character(len=16) :: transform, std_name

        transform = "LEVEL"
        if (present(iv_transform)) transform = uppercase(trim(iv_transform))
        do_standardize_extra = .false.
        if (present(standardize_extra_by_prior_iv)) do_standardize_extra = standardize_extra_by_prior_iv
        if (is_change_transform(transform)) then
            if (present(horizons)) then
                allocate(use_horizons(size(horizons)))
                use_horizons = horizons
            else
                hvals = [1]
                allocate(use_horizons(1))
                use_horizons = hvals
            end if
        else
            allocate(use_horizons(1))
            use_horizons = [0]
        end if

        write(*,'(/,A)') "Implied-volatility correlations:"
        write(*,'(A,A)') "Only mapped implied-vol columns are used from ", trim(implied_vol_file)
        write(*,'(A,A)') "IV transform: ", trim(transform)
        write(*,'(A)') "IV series are lagged one observation: GARCH forecast date D is matched to IV from the prior IV date."
        if (present(extra_names) .and. present(extra_dates) .and. present(extra_values)) then
            write(*,'(A)') "Extra non-forecast rows are matched to same-day IV changes, not lagged IV."
            if (present(extra_signs)) then
                if (any(extra_signs < 0.0_dp)) then
                    write(*,'(A)') "Extra rows with negative sign are sign-flipped before correlation."
                end if
            end if
            if (do_standardize_extra) then
                write(*,'(A)') "STD_* rows divide h-day log returns by prior IV scaled to the h-day horizon."
            end if
        end if
        if (is_change_transform(transform)) then
            write(*,'(A)') "For DIFF/LOG_DIFF, horizon is the h-day change used for both IV and GARCH forecast vol."
        else if (present(horizons)) then
            write(*,'(A)') "Supplied horizons are ignored for LEVEL/LOG transforms."
        end if
        write(*,'(A)') "Asset     IV_index Transform horizon Model              n       corr corr_rank first_date  last_date"
        write(*,'(A)') repeat("-", 100)

        do imap = 1, size(map_assets)
            iasset = findloc(uppercase(asset_names), uppercase(map_assets(imap)), dim=1)
            iiv = findloc(uppercase(iv_names), uppercase(map_iv_indices(imap)), dim=1)
            if (iasset == 0) then
                write(*,'(A,A,A)') "Skipping mapped asset not fit or not found: ", &
                    trim(map_assets(imap)), " (not in selected price columns)"
                cycle
            end if
            if (iiv == 0) then
                write(*,'(A,A)') "Skipping missing implied-vol index: ", trim(map_iv_indices(imap))
                cycle
            end if

            do ih = 1, size(use_horizons)
                if (use_horizons(ih) < 0) cycle
                call transform_series(iv_dates, iv_values(:,iiv), transform, use_horizons(ih), &
                                      iv_corr_dates, iv_corr_values)
                iv_lag_dates = iv_corr_dates
                iv_lag_values = iv_corr_values
                call lag_series_one(iv_lag_dates, iv_lag_values)

                max_rows = count(fit_model_have)
                if (trim(transform) == "LOG_DIFF" .and. present(extra_names) .and. &
                    present(extra_dates) .and. present(extra_values)) then
                    max_rows = max_rows + size(extra_names)
                    if (do_standardize_extra) max_rows = max_rows + size(extra_names)
                end if
                allocate(row_model(max_rows), row_corr(max_rows), row_nmatch(max_rows), &
                         row_first_index(max_rows), row_last_index(max_rows), row_valid(max_rows))
                nrow = 0

                do imodel = 1, size(model_names)
                    if (.not. fit_model_have(imodel)) cycle
                    call transform_series(forecast_dates, vol_forecasts(:,iasset,imodel), transform, &
                                          use_horizons(ih), garch_corr_dates, garch_corr_values)
                    call aligned_corr(garch_corr_dates, garch_corr_values, iv_lag_dates, iv_lag_values, &
                                      corr, nmatch, first_index, last_index)
                    nrow = nrow + 1
                    row_model(nrow) = trim(model_names(imodel))
                    row_corr(nrow) = corr
                    row_nmatch(nrow) = nmatch
                    row_first_index(nrow) = first_index
                    row_last_index(nrow) = last_index
                    row_valid(nrow) = nmatch > 1
                    deallocate(garch_corr_dates, garch_corr_values)
                end do

                if (trim(transform) == "LOG_DIFF" .and. present(extra_names) .and. &
                    present(extra_dates) .and. present(extra_values)) then
                    do iextra = 1, size(extra_names)
                        extra_sign = 1.0_dp
                        if (present(extra_signs)) extra_sign = extra_signs(iextra)
                        call transform_series(extra_dates, extra_values(:,iasset,iextra), transform, &
                                              use_horizons(ih), garch_corr_dates, garch_corr_values)
                        garch_corr_values = extra_sign * garch_corr_values
                        call aligned_corr(garch_corr_dates, garch_corr_values, iv_corr_dates, iv_corr_values, &
                                          corr, nmatch, first_index, last_index)
                        nrow = nrow + 1
                        row_model(nrow) = trim(extra_names(iextra))
                        row_corr(nrow) = corr
                        row_nmatch(nrow) = nmatch
                        row_first_index(nrow) = first_index
                        row_last_index(nrow) = last_index
                        row_valid(nrow) = nmatch > 1
                        deallocate(garch_corr_dates, garch_corr_values)

                        if (do_standardize_extra) then
                            call standardized_log_return_by_prior_iv(extra_dates, extra_values(:,iasset,iextra), &
                                                                     iv_dates, iv_values(:,iiv), use_horizons(ih), &
                                                                     std_dates, std_values)
                            std_values = extra_sign * std_values
                            call aligned_corr(std_dates, std_values, iv_corr_dates, iv_corr_values, &
                                              corr, nmatch, first_index, last_index)
                            std_name = "STD_" // trim(extra_names(iextra))
                            nrow = nrow + 1
                            row_model(nrow) = trim(std_name)
                            row_corr(nrow) = corr
                            row_nmatch(nrow) = nmatch
                            row_first_index(nrow) = first_index
                            row_last_index(nrow) = last_index
                            row_valid(nrow) = nmatch > 1
                            deallocate(std_dates, std_values)
                        end if
                    end do
                end if
                do irow = 1, nrow
                    if (row_valid(irow)) then
                        write(*,'(A9,1X,A8,1X,A9,I8,1X,A16,I7,F11.4,I10,2(1X,A10))') trim(asset_names(iasset)), &
                            trim(iv_names(iiv)), trim(transform), use_horizons(ih), trim(row_model(irow)), &
                            row_nmatch(irow), row_corr(irow), corr_rank(row_corr(irow), row_corr(1:nrow), row_valid(1:nrow)), &
                            date_label(row_first_index(irow)), date_label(row_last_index(irow))
                    else
                        write(*,'(A9,1X,A8,1X,A9,I8,1X,A16,A)') trim(asset_names(iasset)), trim(iv_names(iiv)), &
                            trim(transform), use_horizons(ih), trim(row_model(irow)), &
                            " insufficient overlapping dates"
                    end if
                end do
                deallocate(row_model, row_corr, row_nmatch, row_first_index, row_last_index, row_valid)
                deallocate(iv_corr_dates, iv_corr_values, iv_lag_dates, iv_lag_values)
            end do
        end do
        deallocate(use_horizons)
    end subroutine print_implied_vol_correlations

    integer function corr_rank(value, values, valid)
        real(dp), intent(in) :: value, values(:)
        logical, intent(in) :: valid(:)
        integer :: i

        corr_rank = 1
        do i = 1, size(values)
            if (valid(i) .and. values(i) > value) corr_rank = corr_rank + 1
        end do
    end function corr_rank

    subroutine standardized_log_return_by_prior_iv(price_dates, prices, iv_dates, iv_levels, horizon, out_dates, out_values)
        integer, intent(in) :: price_dates(:), iv_dates(:), horizon
        real(dp), intent(in) :: prices(:), iv_levels(:)
        integer, allocatable, intent(out) :: out_dates(:)
        real(dp), allocatable, intent(out) :: out_values(:)
        integer :: h, n, t, niv, count, iv_idx
        integer, allocatable :: tmp_dates(:)
        real(dp), allocatable :: tmp_values(:)
        real(dp) :: denom

        h = max(horizon, 1)
        n = size(prices)
        niv = size(iv_levels)
        allocate(tmp_dates(max(n-h, 0)), tmp_values(max(n-h, 0)))
        count = 0
        iv_idx = 1
        do t = h + 1, n
            do while (iv_idx <= niv)
                if (iv_dates(iv_idx) >= price_dates(t-h)) exit
                iv_idx = iv_idx + 1
            end do
            if (iv_idx > niv) exit
            if (iv_dates(iv_idx) /= price_dates(t-h)) cycle
            denom = (iv_levels(iv_idx) / 100.0_dp) * sqrt(real(h, dp) / 252.0_dp)
            if (denom <= 0.0_dp) cycle
            count = count + 1
            tmp_dates(count) = price_dates(t)
            tmp_values(count) = (log(prices(t)) - log(prices(t-h))) / denom
        end do

        allocate(out_dates(count), out_values(count))
        out_dates = tmp_dates(1:count)
        out_values = tmp_values(1:count)
        deallocate(tmp_dates, tmp_values)
    end subroutine standardized_log_return_by_prior_iv

end module vol_forecast_compare_mod

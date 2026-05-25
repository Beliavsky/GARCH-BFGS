! Fit MCS-GARCH intraday volatility models to multiple OHLCV price files.
!
! Command-line arguments are interpreted as input CSV files. Arguments containing
! '*' or '?' are expanded with `cmd /c dir /b /s`, so Windows-style globs such as
! c:\python\intraday_prices\*.csv can be used.

module fit_mcsgarch_intraday_batch_mod
    use kind_mod, only: dp
    use date_mod, only: yyyymmdd
    use market_data_mod, only: ohlcv_series_t, read_intraday_prices_csv, filter_intraday_session, intraday_bin_ids
    use garch_mcsgarch_mod, only: mcsgarch_fit_result_t, fit_mcsgarch, fit_mcsgarch_nagarch, fit_mcsgarch_gjr, &
                                  fit_mcsgarch_t, fit_mcsgarch_nagarch_t, fit_mcsgarch_gjr_t, &
                                  fit_mcsgarch_fs_skewt, fit_mcsgarch_nagarch_fs_skewt, fit_mcsgarch_gjr_fs_skewt
    use stats_mod, only: mean
    implicit none
    private

    integer, parameter :: max_files = 256
    integer, parameter :: path_len = 512
    real(dp), parameter :: min_daily_var = 1.0e-12_dp
    logical, parameter :: smooth_diurnal_curve = .true.
    integer, parameter :: diurnal_smooth_half_width = 2
    character(len=*), parameter :: default_file_pattern = "c:\python\intraday_prices\*.csv"
    character(len=12), parameter :: model_names(*) = [character(len=12) :: &
        "MCSGARCH", "MCSNAGARCH", "MCSGJRGARCH"]
    character(len=8), parameter :: dist_names(*) = [character(len=8) :: &
        "NORMAL", "T", "FS_SKEWT"]

    type :: batch_result_t
        character(len=path_len) :: filename = ""
        integer :: n_bars = 0
        integer :: n_returns = 0
        integer :: n_bins = 0
        integer :: n_converged = 0
        character(len=36) :: best_aic_model = ""
        character(len=8) :: best_aic_dist = ""
        character(len=36) :: best_bic_model = ""
        character(len=8) :: best_bic_dist = ""
        real(dp) :: best_aic_loglik = 0.0_dp
        real(dp) :: best_bic_loglik = 0.0_dp
        real(dp) :: best_aic = 0.0_dp
        real(dp) :: best_bic = 0.0_dp
        real(dp) :: read_sec = 0.0_dp
        real(dp) :: fit_sec = 0.0_dp
        real(dp) :: elapsed_sec = 0.0_dp
        logical :: ok = .false.
    end type batch_result_t

    public :: run_fit_mcsgarch_intraday_batch

contains

    ! Fit all configured MCS-GARCH models to each requested intraday price file.
    subroutine run_fit_mcsgarch_intraday_batch()
        character(len=path_len), allocatable :: files(:)
        type(batch_result_t), allocatable :: results(:)
        integer :: nfiles, i
        real(dp) :: t0, t1

        call collect_input_files(files, nfiles)
        if (nfiles < 1) error stop "xfit_mcsgarch_intraday_batch: no input files found"
        allocate(results(nfiles))

        print '(A)', "MCS-GARCH intraday batch fits"
        print '(A,I0)', "Input files: ", nfiles
        print '(A,I0,A,I0)', "Configured models: ", size(model_names), "  distributions: ", size(dist_names)
        print '(A)', ""

        call cpu_time(t0)
        do i = 1, nfiles
            call fit_one_file(trim(files(i)), results(i))
        end do
        call cpu_time(t1)

        call print_batch_summary(results)
        print '(A,F10.3)', "Batch elapsed seconds: ", t1 - t0
        deallocate(files, results)
    end subroutine run_fit_mcsgarch_intraday_batch

    ! Expand command-line filenames and wildcard patterns into a file list.
    subroutine collect_input_files(files, nfiles)
        character(len=path_len), allocatable, intent(out) :: files(:)
        integer, intent(out) :: nfiles
        character(len=path_len) :: arg
        integer :: nargs, i

        allocate(files(max_files))
        nfiles = 0
        nargs = command_argument_count()
        if (nargs < 1) then
            call append_pattern_files(default_file_pattern, files, nfiles)
        else
            do i = 1, nargs
                call get_command_argument(i, arg)
                if (has_wildcard(trim(arg))) then
                    call append_pattern_files(trim(arg), files, nfiles)
                else
                    call append_file(trim(arg), files, nfiles)
                end if
            end do
        end if
        files = files(1:nfiles)
    end subroutine collect_input_files

    ! Return true when a command argument contains wildcard metacharacters.
    logical function has_wildcard(text)
        character(len=*), intent(in) :: text

        has_wildcard = index(text, "*") > 0 .or. index(text, "?") > 0
    end function has_wildcard

    ! Append one filename to the fixed-capacity batch list.
    subroutine append_file(filename, files, nfiles)
        character(len=*), intent(in) :: filename
        character(len=path_len), intent(inout) :: files(:)
        integer, intent(inout) :: nfiles

        if (len_trim(filename) < 1) return
        if (nfiles >= size(files)) error stop "append_file: too many input files"
        nfiles = nfiles + 1
        files(nfiles) = filename
    end subroutine append_file

    ! Expand a Windows wildcard pattern and append matching paths.
    subroutine append_pattern_files(pattern, files, nfiles)
        character(len=*), intent(in) :: pattern
        character(len=path_len), intent(inout) :: files(:)
        integer, intent(inout) :: nfiles
        character(len=*), parameter :: tmp_file = "xfit_mcsgarch_intraday_batch_files.tmp"
        character(len=3*path_len) :: command
        character(len=path_len) :: line
        integer :: unit, ios, exitstat

        open(newunit=unit, file=tmp_file, status="replace", action="write")
        close(unit)
        command = 'cmd /c dir /b /s "' // trim(pattern) // '" > "' // tmp_file // '"'
        call execute_command_line(trim(command), exitstat=exitstat)
        if (exitstat /= 0) return

        open(newunit=unit, file=tmp_file, status="old", action="read")
        do
            read(unit,'(A)',iostat=ios) line
            if (ios /= 0) exit
            call append_file(trim(line), files, nfiles)
        end do
        close(unit, status="delete")
    end subroutine append_pattern_files

    ! Read one file, build inputs, fit all model/distribution pairs, and summarize selection.
    subroutine fit_one_file(filename, result)
        character(len=*), intent(in) :: filename
        type(batch_result_t), intent(out) :: result
        type(ohlcv_series_t) :: bars, regular_bars
        type(mcsgarch_fit_result_t), allocatable :: fit(:)
        real(dp), allocatable :: returns(:), daily_var(:), diurnal_fit(:,:), q_fit(:,:), fit_sec(:)
        integer, allocatable :: bin_id(:), return_dates(:)
        integer :: max_iter, nfit, ifit, imodel, idist
        real(dp) :: gtol, t_start, t0, t1, t_end

        result%filename = filename
        max_iter = 120
        gtol = 1.0e-5_dp
        call cpu_time(t_start)

        print '(A,A)', "Fitting: ", trim(filename)
        call cpu_time(t0)
        call read_intraday_prices_csv(trim(filename), bars)
        call cpu_time(t1)
        result%read_sec = t1 - t0

        call filter_intraday_session(bars, regular_bars)
        call build_mcsgarch_inputs(regular_bars, returns, daily_var, bin_id, return_dates)
        result%n_bars = regular_bars%nobs()
        result%n_returns = size(returns)
        result%n_bins = count_populated_bins(bin_id)

        nfit = size(model_names) * size(dist_names)
        allocate(fit(nfit), fit_sec(nfit), diurnal_fit(size(returns), nfit), q_fit(size(returns), nfit))

        ifit = 0
        do idist = 1, size(dist_names)
            do imodel = 1, size(model_names)
                ifit = ifit + 1
                call cpu_time(t0)
                call fit_configured_model(trim(model_names(imodel)), trim(dist_names(idist)), &
                                          returns, daily_var, bin_id, max_iter, gtol, &
                                          fit(ifit), diurnal_fit(:, ifit), q_fit(:, ifit))
                call cpu_time(t1)
                fit_sec(ifit) = t1 - t0
            end do
        end do

        call summarize_fit_selection(returns, bin_id, fit, fit_sec, result)
        call cpu_time(t_end)
        result%elapsed_sec = t_end - t_start
        result%ok = .true.
        print '(A,A,A,A,A,A,A,A,A,F10.3)', "  AIC: ", trim(result%best_aic_model), " ", &
              trim(result%best_aic_dist), "  BIC: ", trim(result%best_bic_model), " ", &
              trim(result%best_bic_dist), "  elapsed: ", result%elapsed_sec

        deallocate(returns, daily_var, bin_id, return_dates, fit, fit_sec, diurnal_fit, q_fit)
    end subroutine fit_one_file

    ! Fit one configured dynamic model and innovation distribution.
    subroutine fit_configured_model(model_name, dist_name, returns, daily_var, bin_id, max_iter, gtol, fit, diurnal_var, q)
        character(len=*), intent(in) :: model_name, dist_name
        real(dp), intent(in) :: returns(:), daily_var(:), gtol
        integer, intent(in) :: bin_id(:), max_iter
        type(mcsgarch_fit_result_t), intent(out) :: fit
        real(dp), intent(out) :: diurnal_var(:), q(:)

        select case (trim(model_name))
        case ("MCSGARCH")
            if (trim(dist_name) == "T") then
                call fit_mcsgarch_t(returns, daily_var, bin_id, max_iter, gtol, fit, diurnal_var, q, &
                                    smooth_diurnal=smooth_diurnal_curve, &
                                    smooth_half_width=diurnal_smooth_half_width)
            else if (trim(dist_name) == "FS_SKEWT") then
                call fit_mcsgarch_fs_skewt(returns, daily_var, bin_id, max_iter, gtol, fit, diurnal_var, q, &
                                           smooth_diurnal=smooth_diurnal_curve, &
                                           smooth_half_width=diurnal_smooth_half_width)
            else
                call fit_mcsgarch(returns, daily_var, bin_id, max_iter, gtol, fit, diurnal_var, q, &
                                  smooth_diurnal=smooth_diurnal_curve, &
                                  smooth_half_width=diurnal_smooth_half_width)
            end if
        case ("MCSNAGARCH")
            if (trim(dist_name) == "T") then
                call fit_mcsgarch_nagarch_t(returns, daily_var, bin_id, max_iter, gtol, fit, diurnal_var, q, &
                                            smooth_diurnal=smooth_diurnal_curve, &
                                            smooth_half_width=diurnal_smooth_half_width)
            else if (trim(dist_name) == "FS_SKEWT") then
                call fit_mcsgarch_nagarch_fs_skewt(returns, daily_var, bin_id, max_iter, gtol, fit, diurnal_var, q, &
                                                   smooth_diurnal=smooth_diurnal_curve, &
                                                   smooth_half_width=diurnal_smooth_half_width)
            else
                call fit_mcsgarch_nagarch(returns, daily_var, bin_id, max_iter, gtol, fit, diurnal_var, q, &
                                          smooth_diurnal=smooth_diurnal_curve, &
                                          smooth_half_width=diurnal_smooth_half_width)
            end if
        case ("MCSGJRGARCH")
            if (trim(dist_name) == "T") then
                call fit_mcsgarch_gjr_t(returns, daily_var, bin_id, max_iter, gtol, fit, diurnal_var, q, &
                                        smooth_diurnal=smooth_diurnal_curve, &
                                        smooth_half_width=diurnal_smooth_half_width)
            else if (trim(dist_name) == "FS_SKEWT") then
                call fit_mcsgarch_gjr_fs_skewt(returns, daily_var, bin_id, max_iter, gtol, fit, diurnal_var, q, &
                                               smooth_diurnal=smooth_diurnal_curve, &
                                               smooth_half_width=diurnal_smooth_half_width)
            else
                call fit_mcsgarch_gjr(returns, daily_var, bin_id, max_iter, gtol, fit, diurnal_var, q, &
                                      smooth_diurnal=smooth_diurnal_curve, &
                                      smooth_half_width=diurnal_smooth_half_width)
            end if
        case default
            error stop "fit_configured_model: unsupported model"
        end select
    end subroutine fit_configured_model

    ! Convert regular-session bars into within-session returns and daily variance proxies.
    subroutine build_mcsgarch_inputs(bars, returns, daily_var, bin_id, return_dates)
        type(ohlcv_series_t), intent(in) :: bars
        real(dp), allocatable, intent(out) :: returns(:), daily_var(:)
        integer, allocatable, intent(out) :: bin_id(:), return_dates(:)
        integer, allocatable :: bar_bins(:), day_index(:), day_dates(:), day_first(:), day_last(:)
        real(dp), allocatable :: intraday_rv(:), overnight_sq(:), daily_proxy(:)
        real(dp), allocatable :: returns_all(:), daily_var_all(:)
        integer, allocatable :: bin_all(:), dates_all(:)
        integer :: n, ndays, i, k, d
        real(dp) :: fallback_daily_var

        n = bars%nobs()
        if (n < 3) error stop "build_mcsgarch_inputs: not enough regular-session bars"
        call intraday_bin_ids(bars, bar_bins)
        call map_days(bars, day_index, day_dates, day_first, day_last)
        ndays = size(day_dates)
        allocate(intraday_rv(ndays), overnight_sq(ndays), daily_proxy(ndays))
        intraday_rv = 0.0_dp
        overnight_sq = 0.0_dp

        do i = 2, n
            if (day_index(i) == day_index(i - 1)) then
                intraday_rv(day_index(i)) = intraday_rv(day_index(i)) + log(bars%close(i) / bars%close(i - 1))**2
            end if
        end do
        do d = 2, ndays
            overnight_sq(d) = log(bars%open(day_first(d)) / bars%close(day_last(d - 1)))**2
        end do

        fallback_daily_var = max(mean(intraday_rv(max(1, min(2, ndays)):ndays)), min_daily_var)
        daily_proxy(1) = fallback_daily_var
        do d = 2, ndays
            daily_proxy(d) = max(intraday_rv(d - 1) + overnight_sq(d), min_daily_var)
        end do

        allocate(returns_all(n - 1), daily_var_all(n - 1), bin_all(n - 1), dates_all(n - 1))
        k = 0
        do i = 2, n
            if (day_index(i) /= day_index(i - 1)) cycle
            k = k + 1
            returns_all(k) = log(bars%close(i) / bars%close(i - 1))
            daily_var_all(k) = daily_proxy(day_index(i))
            bin_all(k) = bar_bins(i)
            dates_all(k) = day_dates(day_index(i))
        end do
        if (k < 3) error stop "build_mcsgarch_inputs: not enough intraday returns"
        allocate(returns(k), daily_var(k), bin_id(k), return_dates(k))
        returns = returns_all(1:k)
        daily_var = daily_var_all(1:k)
        bin_id = bin_all(1:k)
        return_dates = dates_all(1:k)

        deallocate(bar_bins, day_index, day_dates, day_first, day_last, intraday_rv, overnight_sq, daily_proxy)
        deallocate(returns_all, daily_var_all, bin_all, dates_all)
    end subroutine build_mcsgarch_inputs

    ! Map each bar to a trading day and record the first and last bar of each day.
    subroutine map_days(bars, day_index, day_dates, day_first, day_last)
        type(ohlcv_series_t), intent(in) :: bars
        integer, allocatable, intent(out) :: day_index(:), day_dates(:), day_first(:), day_last(:)
        integer :: n, i, ndays, current_date

        n = bars%nobs()
        allocate(day_index(n), day_dates(n), day_first(n), day_last(n))
        ndays = 0
        current_date = -1
        do i = 1, n
            if (yyyymmdd(bars%timestamp(i)%date) /= current_date) then
                ndays = ndays + 1
                current_date = yyyymmdd(bars%timestamp(i)%date)
                day_dates(ndays) = current_date
                day_first(ndays) = i
                if (ndays > 1) day_last(ndays - 1) = i - 1
            end if
            day_index(i) = ndays
        end do
        day_last(ndays) = n
        day_dates = day_dates(1:ndays)
        day_first = day_first(1:ndays)
        day_last = day_last(1:ndays)
    end subroutine map_days

    ! Pick the AIC and BIC winners among fitted model/distribution combinations.
    subroutine summarize_fit_selection(returns, bin_id, fit, fit_sec, result)
        real(dp), intent(in) :: returns(:), fit_sec(:)
        integer, intent(in) :: bin_id(:)
        type(mcsgarch_fit_result_t), intent(in) :: fit(:)
        type(batch_result_t), intent(inout) :: result
        real(dp) :: aic, bic, best_aic, best_bic
        integer :: ifit, imodel, idist, k, nobs
        logical :: first

        nobs = size(returns)
        result%n_converged = count(fit%converged)
        result%fit_sec = sum(fit_sec)
        first = .true.
        ifit = 0
        do idist = 1, size(dist_names)
            do imodel = 1, size(model_names)
                ifit = ifit + 1
                k = count_populated_bins(bin_id) + model_param_count(trim(model_names(imodel))) + &
                    dist_param_count(trim(dist_names(idist)))
                aic = -2.0_dp*fit(ifit)%loglik + 2.0_dp*real(k, dp)
                bic = -2.0_dp*fit(ifit)%loglik + real(k, dp)*log(real(nobs, dp))
                if (first .or. aic < best_aic) then
                    best_aic = aic
                    result%best_aic = aic
                    result%best_aic_loglik = fit(ifit)%loglik
                    result%best_aic_model = model_table_label(trim(model_names(imodel)))
                    result%best_aic_dist = trim(dist_names(idist))
                end if
                if (first .or. bic < best_bic) then
                    best_bic = bic
                    result%best_bic = bic
                    result%best_bic_loglik = fit(ifit)%loglik
                    result%best_bic_model = model_table_label(trim(model_names(imodel)))
                    result%best_bic_dist = trim(dist_names(idist))
                end if
                first = .false.
            end do
        end do
    end subroutine summarize_fit_selection

    ! Return the dynamic-parameter count for an MCS model name.
    integer function model_param_count(model_name)
        character(len=*), intent(in) :: model_name

        select case (trim(model_name))
        case ("MCSGARCH")
            model_param_count = 3
        case ("MCSNAGARCH", "MCSGJRGARCH")
            model_param_count = 4
        case default
            error stop "model_param_count: unsupported model"
        end select
    end function model_param_count

    ! Return the innovation-distribution parameter count.
    integer function dist_param_count(dist_name)
        character(len=*), intent(in) :: dist_name

        select case (trim(dist_name))
        case ("NORMAL")
            dist_param_count = 0
        case ("T")
            dist_param_count = 1
        case ("FS_SKEWT")
            dist_param_count = 2
        case default
            error stop "dist_param_count: unsupported distribution"
        end select
    end function dist_param_count

    ! Return display label used in the summary table.
    character(len=36) function model_table_label(model_name)
        character(len=*), intent(in) :: model_name

        select case (trim(model_name))
        case ("MCSGARCH")
            model_table_label = "MCS-GARCH"
        case ("MCSNAGARCH")
            model_table_label = "MCS-NAGARCH"
        case ("MCSGJRGARCH")
            model_table_label = "MCS-GJRGARCH"
        case default
            model_table_label = trim(model_name)
        end select
    end function model_table_label

    ! Count the number of intraday bins represented in the return sample.
    integer function count_populated_bins(bin_id)
        integer, intent(in) :: bin_id(:)
        integer, allocatable :: bin_count(:)
        integer :: i

        allocate(bin_count(maxval(bin_id)))
        bin_count = 0
        do i = 1, size(bin_id)
            bin_count(bin_id(i)) = bin_count(bin_id(i)) + 1
        end do
        count_populated_bins = count(bin_count > 0)
        deallocate(bin_count)
    end function count_populated_bins

    ! Print a compact summary across fitted input files.
    subroutine print_batch_summary(results)
        type(batch_result_t), intent(in) :: results(:)
        integer :: i

        print '(A)', ""
        print '(A)', "MCS-GARCH intraday batch summary"
        print '(A)', "--------------------------------------------------------------------------------------------------------------------------------------------------------------"
        print '(A32,1X,A8,1X,A8,1X,A5,1X,A6,1X,A18,1X,A8,1X,A18,1X,A8,1X,A12,1X,A12,1X,A8,1X,A8,1X,A8)', &
              "File", "n_bars", "n_ret", "bins", "conv", "AIC_model", "AIC_dist", "BIC_model", "BIC_dist", &
              "AIC", "BIC", "read_s", "fit_s", "total_s"
        print '(A)', "--------------------------------------------------------------------------------------------------------------------------------------------------------------"
        do i = 1, size(results)
            call print_batch_summary_row(results(i))
        end do
        print '(A)', "--------------------------------------------------------------------------------------------------------------------------------------------------------------"
    end subroutine print_batch_summary

    ! Print one batch summary row.
    subroutine print_batch_summary_row(result)
        type(batch_result_t), intent(in) :: result

        print '(A32,1X,I8,1X,I8,1X,I5,1X,I6,1X,A18,1X,A8,1X,A18,1X,A8,1X,F12.3,1X,F12.3,1X,F8.3,1X,F8.3,1X,F8.3)', &
              short_filename(result%filename), result%n_bars, result%n_returns, result%n_bins, result%n_converged, &
              result%best_aic_model, result%best_aic_dist, result%best_bic_model, result%best_bic_dist, &
              result%best_aic, result%best_bic, result%read_sec, result%fit_sec, result%elapsed_sec
    end subroutine print_batch_summary_row

    ! Return the final path component for compact table display.
    function short_filename(path) result(name)
        character(len=*), intent(in) :: path
        character(len=32) :: name
        integer :: pos1, pos2, pos

        pos1 = scan(trim(path), "\", back=.true.)
        pos2 = scan(trim(path), "/", back=.true.)
        pos = max(pos1, pos2)
        if (pos > 0) then
            name = path(pos + 1:)
        else
            name = path
        end if
    end function short_filename

end module fit_mcsgarch_intraday_batch_mod

program xfit_mcsgarch_intraday_batch
    use fit_mcsgarch_intraday_batch_mod, only: run_fit_mcsgarch_intraday_batch
    implicit none
    call run_fit_mcsgarch_intraday_batch()
end program xfit_mcsgarch_intraday_batch

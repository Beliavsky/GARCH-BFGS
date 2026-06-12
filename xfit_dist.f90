! Fit iid distributions to numeric series read from a CSV file.
!
! The CSV may have a header row and may use the first column as an integer row
! index.  Index limits are inclusive and make it easy to fit subsets such as
! row ranges or YYYYMMDD date ranges.  Columns after the optional index column
! are treated as iid numeric series.
!
! Usage: xfit_dist.exe [input_csv] [max_rows] [max_cols] [max_iter=N] [stats=true|false] [init=fixed|stats|robust|warm]
! Use max_rows <= 0 or max_cols <= 0 for no limit.

program xfit_dist
    use date_mod, only: print_program_header
    use kind_mod, only: dp
    use csv_mod, only: read_numeric_csv
    use stats_mod, only: mean, sd, column_summary_stats, column_robust_tail_ratios, print_column_summary_stats
    use distributions_mod, only: dist_id_from_name, dist_fixed_shape_from_name, dist_npar_std, fit_dist, fit_dist_std, &
                                 dist_normal, dist_t, dist_fs_skewt, dist_warm_shape_start, &
                                 dist_warm_shape_start_scaled, dist_warm_sigma_start
    use strings_mod, only: uppercase
    use rank_mod, only: rank_asc
    implicit none

    character(len=*), parameter :: default_file = "dist_input.csv"
    logical, parameter :: has_header = .true.
    logical, parameter :: has_index_col = .true.
    logical, parameter :: use_index_min = .false.
    logical, parameter :: use_index_max = .false.
    logical, parameter :: standardized_input = .false.
    integer, parameter :: index_min = 0
    integer, parameter :: index_max = huge(0)
    integer, parameter :: default_max_rows = 0
    integer, parameter :: default_max_cols = 1000
    integer, parameter :: default_max_iter = 500
    character(len=16), parameter :: dist_list(*) = [character(len=16) :: &
        "NORMAL", "LOGISTIC", "LAPLACE", "SECH", "T_4", "T_6", "T_8", "T", "GED", "NIG_SYM", "NIG", "FS_SKEWT"]

    character(len=256) :: input_file, arg, opt_name, opt_value
    integer, allocatable :: row_index(:)
    integer, allocatable :: row_dist_id(:), row_iter(:), aic_rank(:), bic_rank(:)
    character(len=16), allocatable :: row_dist_name(:)
    character(len=32), allocatable :: col_names(:)
    character(len=12), allocatable :: best_aic_1(:), best_aic_2(:), best_bic_1(:), best_bic_2(:)
    real(dp), allocatable :: values(:,:), x(:)
    real(dp), allocatable :: row_mu(:), row_sigma(:), row_shape(:), row_xi(:), row_loglik(:), row_aic(:), row_bic(:)
    real(dp), allocatable :: row_sec(:)
    real(dp), allocatable :: best_aic_diff(:), best_bic_diff(:)
    real(dp), allocatable :: stat_median(:), stat_mean(:), stat_sd(:), stat_skew(:), stat_exkurt(:)
    real(dp), allocatable :: stat_min(:), stat_max(:), stat_tail_ratio(:)
    real(dp) :: mu, sigma, shape, xi, loglik, aic, bic, fit_start, fit_end, fit_sec, t_start, t_end
    real(dp) :: best_simple_aic, warm_shape, warm_sigma, fixed_shape
    real(dp) :: read_start, read_end, stats_timer_start, stats_timer_end, read_sec, stats_sec, fit_total_sec
    integer :: nobs, ncols, icol, idist, ifit, nfit, dist_id, nparam, niter, max_rows, max_cols, max_iter
    integer :: best_simple_dist_id
    integer :: best1, best2
    integer :: npos
    integer :: ios, iarg, eq_pos
    logical :: converged, print_stats, stats_start, robust_start, warm_start, has_fixed_shape
    logical, allocatable :: row_converged(:)

    call print_program_header("xfit_dist.f90")
    call cpu_time(t_start)

    input_file = default_file
    max_rows = default_max_rows
    max_cols = default_max_cols
    max_iter = default_max_iter
    print_stats = .false.
    stats_start = .false.
    robust_start = .false.
    warm_start = .false.

    call get_command_argument(1, input_file)
    if (len_trim(input_file) == 0) input_file = default_file
    npos = 0
    do iarg = 2, command_argument_count()
        call get_command_argument(iarg, arg)
        eq_pos = index(arg, "=")
        if (eq_pos == 0) then
            npos = npos + 1
            if (npos == 1) then
                read(arg, *, iostat=ios) max_rows
                if (ios /= 0) then
                    print '(A)', "Usage: xfit_dist.exe [input_csv] [max_rows] [max_cols] [max_iter=N] [stats=true|false] [init=fixed|stats|robust|warm]"
                    error stop "xfit_dist: invalid max_rows"
                end if
            else if (npos == 2) then
                read(arg, *, iostat=ios) max_cols
                if (ios /= 0) then
                    print '(A)', "Usage: xfit_dist.exe [input_csv] [max_rows] [max_cols] [max_iter=N] [stats=true|false] [init=fixed|stats|robust|warm]"
                    error stop "xfit_dist: invalid max_cols"
                end if
            else
                print '(A)', "Usage: xfit_dist.exe [input_csv] [max_rows] [max_cols] [max_iter=N] [stats=true|false] [init=fixed|stats|robust|warm]"
                error stop "xfit_dist: too many positional arguments"
            end if
        else if (eq_pos < 2) then
            print '(A)', "Usage: xfit_dist.exe [input_csv] [max_rows] [max_cols] [max_iter=N] [stats=true|false] [init=fixed|stats|robust|warm]"
            error stop "xfit_dist: expected name=value option"
        else
            opt_name = uppercase(adjustl(arg(:eq_pos-1)))
            opt_value = adjustl(arg(eq_pos+1:))
            select case (trim(opt_name))
            case ("MAX_ITER")
                read(opt_value, *, iostat=ios) max_iter
                if (ios /= 0) then
                    print '(A)', "Usage: xfit_dist.exe [input_csv] [max_rows] [max_cols] [max_iter=N] [stats=true|false] [init=fixed|stats|robust|warm]"
                    error stop "xfit_dist: invalid max_iter"
                end if
            case ("STATS", "PRINT_STATS")
                select case (trim(uppercase(opt_value)))
                case ("TRUE", "T", ".TRUE.", "1", "YES", "Y")
                    print_stats = .true.
                case ("FALSE", "F", ".FALSE.", "0", "NO", "N")
                    print_stats = .false.
                case default
                    print '(A)', "Usage: xfit_dist.exe [input_csv] [max_rows] [max_cols] [max_iter=N] [stats=true|false] [init=fixed|stats|robust|warm]"
                    error stop "xfit_dist: invalid stats option"
                end select
            case ("INIT", "START")
                select case (trim(uppercase(opt_value)))
                case ("FIXED", "DEFAULT")
                    stats_start = .false.
                    robust_start = .false.
                    warm_start = .false.
                case ("STATS", "SUMMARY")
                    stats_start = .true.
                    robust_start = .false.
                    warm_start = .false.
                case ("ROBUST")
                    stats_start = .false.
                    robust_start = .true.
                    warm_start = .false.
                case ("WARM", "SEQUENTIAL")
                    stats_start = .false.
                    robust_start = .false.
                    warm_start = .true.
                case default
                    print '(A)', "Usage: xfit_dist.exe [input_csv] [max_rows] [max_cols] [max_iter=N] [stats=true|false] [init=fixed|stats|robust|warm]"
                    error stop "xfit_dist: invalid init option"
                end select
            case ("STATS_START")
                select case (trim(uppercase(opt_value)))
                case ("TRUE", "T", ".TRUE.", "1", "YES", "Y")
                    stats_start = .true.
                    robust_start = .false.
                    warm_start = .false.
                case ("FALSE", "F", ".FALSE.", "0", "NO", "N")
                    stats_start = .false.
                case default
                    print '(A)', "Usage: xfit_dist.exe [input_csv] [max_rows] [max_cols] [max_iter=N] [stats=true|false] [init=fixed|stats|robust|warm]"
                    error stop "xfit_dist: invalid stats_start option"
                end select
            case default
                print '(A,A)', "Unsupported option: ", trim(opt_name)
                print '(A)', "Usage: xfit_dist.exe [input_csv] [max_rows] [max_cols] [max_iter=N] [stats=true|false] [init=fixed|stats|robust|warm]"
                error stop "xfit_dist: unsupported option"
            end select
        end if
    end do
    if (max_iter < 1) then
        print '(A)', "Usage: xfit_dist.exe [input_csv] [max_rows] [max_cols] [max_iter=N] [stats=true|false] [init=fixed|stats|robust|warm]"
        error stop "xfit_dist: max_iter must be positive"
    end if
    if (max_rows <= 0) max_rows = huge(max_rows)
    if (max_cols <= 0) max_cols = huge(max_cols)

    call cpu_time(read_start)
    if (use_index_min .and. use_index_max) then
        call read_numeric_csv(input_file, row_index, col_names, values, has_header, has_index_col, &
                              index_min, index_max, max_col=max_cols, max_rows=max_rows)
    else if (use_index_min) then
        call read_numeric_csv(input_file, row_index, col_names, values, has_header, has_index_col, &
                              index_min=index_min, max_col=max_cols, max_rows=max_rows)
    else if (use_index_max) then
        call read_numeric_csv(input_file, row_index, col_names, values, has_header, has_index_col, &
                              index_max=index_max, max_col=max_cols, max_rows=max_rows)
    else
        call read_numeric_csv(input_file, row_index, col_names, values, has_header, has_index_col, &
                              max_col=max_cols, max_rows=max_rows)
    end if
    call cpu_time(read_end)
    read_sec = read_end - read_start

    nobs = size(values, 1)
    ncols = size(values, 2)
    stats_sec = 0.0_dp
    if (print_stats .or. stats_start .or. robust_start) then
        call cpu_time(stats_timer_start)
        call column_summary_stats(values, stat_median, stat_mean, stat_sd, stat_skew, stat_exkurt, stat_min, stat_max)
        call cpu_time(stats_timer_end)
        stats_sec = stats_timer_end - stats_timer_start
    end if
    if (robust_start) then
        call cpu_time(stats_timer_start)
        call column_robust_tail_ratios(values, stat_tail_ratio)
        call cpu_time(stats_timer_end)
        stats_sec = stats_sec + stats_timer_end - stats_timer_start
    end if

    print '(A,A)', "Input file: ", trim(input_file)
    print '(A,I0,A,I0,A,I0,A,I0,/)', "Using ", nobs, " rows for ", ncols, &
        " series; index range ", row_index(1), " to ", row_index(nobs)
    print '(A,I0)', "Maximum fit iterations: ", max_iter
    if (warm_start) then
        print '(A,/)', "Initial guess: warm"
    else if (robust_start) then
        print '(A,/)', "Initial guess: robust"
    else if (stats_start) then
        print '(A,/)', "Initial guess: stats"
    else
        print '(A,/)', "Initial guess: fixed"
    end if

    if (print_stats) then
        call print_column_summary_stats("Column summary statistics", col_names, stat_median, stat_mean, stat_sd, &
                                        stat_skew, stat_exkurt, stat_min, stat_max)
    end if

    if (standardized_input) then
        print '(A)', "Distribution fits to standardized input columns"
    else
        print '(A)', "Distribution fits to raw input columns"
    end if
    print '(A)', repeat("-", 132)
    print '(A16,1X,A12,1X,A10,1X,A10,1X,A10,1X,A10,1X,A12,1X,A12,1X,A12,1X,A8,1X,A8,1X,A4,1X,A6,1X,A10)', &
        "Series", "Dist", "mu", "sigma", "shape", "xi", "logL", "AIC", "BIC", "AIC_rank", "BIC_rank", "conv", "iter", "sec"
    print '(A)', repeat("-", 132)

    allocate(x(nobs))
    allocate(row_dist_id(size(dist_list)), row_dist_name(size(dist_list)), row_mu(size(dist_list)), row_sigma(size(dist_list)), &
             row_shape(size(dist_list)), row_xi(size(dist_list)), row_loglik(size(dist_list)), row_aic(size(dist_list)), &
             row_bic(size(dist_list)), row_sec(size(dist_list)), row_iter(size(dist_list)), &
             aic_rank(size(dist_list)), bic_rank(size(dist_list)), row_converged(size(dist_list)))
    allocate(best_aic_1(ncols), best_aic_2(ncols), best_bic_1(ncols), best_bic_2(ncols), &
             best_aic_diff(ncols), best_bic_diff(ncols))
    fit_total_sec = 0.0_dp
    do icol = 1, ncols
        nfit = 0
        best_simple_aic = huge(best_simple_aic)
        best_simple_dist_id = dist_normal
        do idist = 1, size(dist_list)
            dist_id = dist_id_from_name(dist_list(idist))
            has_fixed_shape = dist_fixed_shape_from_name(dist_list(idist), fixed_shape)
            if (dist_id == 0) then
                print '(A,A)', "Unsupported distribution: ", trim(dist_list(idist))
                cycle
            end if
            if (standardized_input) then
                mu = mean(values(:, icol))
                sigma = sd(values(:, icol))
                x = (values(:, icol) - mu) / max(sigma, 1.0e-20_dp)
                call cpu_time(fit_start)
                if (has_fixed_shape) then
                    call fit_dist_std(x, nobs, dist_id, shape, loglik, converged, niter, max_iter, &
                                      fixed_shape_in=fixed_shape, xi_out=xi)
                else if (warm_start .and. dist_npar_std(dist_id) > 0) then
                    if (dist_id == dist_fs_skewt) then
                        warm_shape = dist_warm_shape_start(dist_t, best_simple_dist_id)
                    else
                        warm_shape = dist_warm_shape_start(dist_id, best_simple_dist_id)
                    end if
                    call fit_dist_std(x, nobs, dist_id, shape, loglik, converged, niter, max_iter, &
                                      start_shape_in=warm_shape, xi_out=xi)
                else if (robust_start) then
                    call fit_dist_std(x, nobs, dist_id, shape, loglik, converged, niter, max_iter, &
                                      start_exkurt_in=stat_exkurt(icol), robust_start_in=robust_start, &
                                      start_tail_ratio_in=stat_tail_ratio(icol), xi_out=xi)
                else if (stats_start) then
                    call fit_dist_std(x, nobs, dist_id, shape, loglik, converged, niter, max_iter, &
                                      stats_start, stat_exkurt(icol), xi_out=xi)
                else
                    call fit_dist_std(x, nobs, dist_id, shape, loglik, converged, niter, max_iter, stats_start, xi_out=xi)
                end if
                call cpu_time(fit_end)
                nparam = merge(0, dist_npar_std(dist_id), has_fixed_shape)
            else
                call cpu_time(fit_start)
                if (has_fixed_shape) then
                    call fit_dist(values(:, icol), nobs, dist_id, mu, sigma, shape, loglik, converged, niter, &
                                  max_iter, fixed_shape_in=fixed_shape, xi_out=xi)
                else if (warm_start .and. dist_npar_std(dist_id) > 0) then
                    if (dist_id == dist_fs_skewt) then
                        warm_shape = dist_warm_shape_start_scaled(dist_t, best_simple_dist_id)
                        warm_sigma = dist_warm_sigma_start(dist_t, best_simple_dist_id)
                    else
                        warm_shape = dist_warm_shape_start_scaled(dist_id, best_simple_dist_id)
                        warm_sigma = dist_warm_sigma_start(dist_id, best_simple_dist_id)
                    end if
                    call fit_dist(values(:, icol), nobs, dist_id, mu, sigma, shape, loglik, converged, niter, &
                                  max_iter, start_shape_in=warm_shape, start_sigma_in=warm_sigma, xi_out=xi)
                else if (robust_start) then
                    call fit_dist(values(:, icol), nobs, dist_id, mu, sigma, shape, loglik, converged, niter, &
                                  max_iter, start_exkurt_in=stat_exkurt(icol), robust_start_in=robust_start, &
                                  start_tail_ratio_in=stat_tail_ratio(icol), xi_out=xi)
                else if (stats_start) then
                    call fit_dist(values(:, icol), nobs, dist_id, mu, sigma, shape, loglik, converged, niter, &
                                  max_iter, stats_start, stat_exkurt(icol), xi_out=xi)
                else
                    call fit_dist(values(:, icol), nobs, dist_id, mu, sigma, shape, loglik, converged, niter, &
                                  max_iter, stats_start, xi_out=xi)
                end if
                call cpu_time(fit_end)
                nparam = 2 + merge(0, dist_npar_std(dist_id), has_fixed_shape)
            end if
            fit_sec = fit_end - fit_start
            fit_total_sec = fit_total_sec + fit_sec
            aic = 2.0_dp*real(nparam, dp) - 2.0_dp*loglik
            bic = log(real(nobs, dp))*real(nparam, dp) - 2.0_dp*loglik
            if (.not. has_fixed_shape .and. dist_npar_std(dist_id) == 0 .and. aic < best_simple_aic) then
                best_simple_aic = aic
                best_simple_dist_id = dist_id
            end if

            nfit = nfit + 1
            row_dist_id(nfit) = dist_id
            row_dist_name(nfit) = dist_list(idist)
            row_mu(nfit) = mu
            row_sigma(nfit) = sigma
            row_shape(nfit) = shape
            row_xi(nfit) = xi
            row_loglik(nfit) = loglik
            row_aic(nfit) = aic
            row_bic(nfit) = bic
            row_sec(nfit) = fit_sec
            row_iter(nfit) = niter
            row_converged(nfit) = converged
        end do
        call rank_asc(row_aic(:nfit), aic_rank(:nfit))
        call rank_asc(row_bic(:nfit), bic_rank(:nfit))
        best1 = 1
        best2 = 2
        if (row_aic(best2) < row_aic(best1)) then
            best1 = 2
            best2 = 1
        end if
        do ifit = 3, nfit
            if (row_aic(ifit) < row_aic(best1)) then
                best2 = best1
                best1 = ifit
            else if (row_aic(ifit) < row_aic(best2)) then
                best2 = ifit
            end if
        end do
        best_aic_1(icol) = trim(row_dist_name(best1))
        best_aic_2(icol) = trim(row_dist_name(best2))
        best_aic_diff(icol) = row_aic(best2) - row_aic(best1)

        best1 = 1
        best2 = 2
        if (row_bic(best2) < row_bic(best1)) then
            best1 = 2
            best2 = 1
        end if
        do ifit = 3, nfit
            if (row_bic(ifit) < row_bic(best1)) then
                best2 = best1
                best1 = ifit
            else if (row_bic(ifit) < row_bic(best2)) then
                best2 = ifit
            end if
        end do
        best_bic_1(icol) = trim(row_dist_name(best1))
        best_bic_2(icol) = trim(row_dist_name(best2))
        best_bic_diff(icol) = row_bic(best2) - row_bic(best1)

        do ifit = 1, nfit
            print '(A16,1X,A12,1X,ES10.3,1X,ES10.3,1X,F10.3,1X,F10.3,1X,F12.3,1X,F12.3,1X,F12.3,1X,I8,1X,I8,1X,L4,1X,I6,1X,F10.3)', &
                trim(col_names(icol)), trim(row_dist_name(ifit)), row_mu(ifit), row_sigma(ifit), &
                row_shape(ifit), row_xi(ifit), row_loglik(ifit), row_aic(ifit), row_bic(ifit), aic_rank(ifit), &
                bic_rank(ifit), row_converged(ifit), row_iter(ifit), row_sec(ifit)
        end do
    end do
    print '(A)', repeat("-", 132)

    print '(/,A)', "Top distributions by information criterion"
    print '(A)', repeat("-", 96)
    print '(A16,1X,A12,1X,A12,1X,A12,1X,A12,1X,A12,1X,A12)', &
        "Series", "AIC_1", "AIC_2", "dAIC", "BIC_1", "BIC_2", "dBIC"
    print '(A)', repeat("-", 96)
    do icol = 1, ncols
        print '(A16,1X,A12,1X,A12,1X,F12.3,1X,A12,1X,A12,1X,F12.3)', &
            trim(col_names(icol)), trim(best_aic_1(icol)), trim(best_aic_2(icol)), best_aic_diff(icol), &
            trim(best_bic_1(icol)), trim(best_bic_2(icol)), best_bic_diff(icol)
    end do
    print '(A)', repeat("-", 96)

    call cpu_time(t_end)
    print '(/,A)', "times in seconds"
    print '(A)', repeat("-", 34)
    print '(A12,1X,A10,1X,A10)', "step", "seconds", "frac"
    print '(A)', repeat("-", 34)
    print '(A12,1X,F10.3,1X,F10.3)', "read", read_sec, read_sec / max(t_end - t_start, 1.0e-20_dp)
    print '(A12,1X,F10.3,1X,F10.3)', "stats", stats_sec, stats_sec / max(t_end - t_start, 1.0e-20_dp)
    print '(A12,1X,F10.3,1X,F10.3)', "fit", fit_total_sec, fit_total_sec / max(t_end - t_start, 1.0e-20_dp)
    print '(A12,1X,F10.3,1X,F10.3)', "elapsed", t_end - t_start, 1.0_dp
    print '(A)', repeat("-", 34)
end program xfit_dist

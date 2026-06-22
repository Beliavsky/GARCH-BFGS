! Sample statistics for real(dp) arrays.
!
! mean: sample mean
! variance: sample variance (denominator n-1)
! sd:       sample standard deviation (denominator n-1)
! sort_real: ascending in-place sort
! covariance: sample covariance of two vectors
! correlation: Pearson correlation of two vectors
! demean_first: subtract mean of first n values from a vector
! autocorr:  autocorrelations at lags 1:size(acf)
! print_acf_table: formatted table of ACF rows
! column_summary_stats: compute summary statistics for matrix columns
! column_robust_tail_ratios: quantile tail ratio by matrix column
! print_column_summary_stats: formatted summary statistics for matrix columns
! nonnegative_least_squares: coordinate-descent NNLS for y ~= X*b, b >= 0

module stats_mod
    use kind_mod, only: dp
    implicit none
    private
    public :: mean, variance, sd, median, skewness, excess_kurtosis, covariance, correlation, demean_first
    public :: sort_real, autocorr, print_acf_table, column_summary_stats, column_robust_tail_ratios
    public :: print_column_summary_stats, covariance_matrix, correlation_matrix, nonnegative_least_squares

contains

    pure function mean(x) result(m)
        ! Sample mean of x.
        real(dp), intent(in) :: x(:)  ! input array
        real(dp) :: m
        m = sum(x) / size(x)
    end function mean

    pure function variance(x) result(v)
        ! Sample variance of x (denominator n-1).
        real(dp), intent(in) :: x(:)  ! input array
        real(dp) :: v
        integer :: n

        n = size(x)
        if (n > 1) then
            v = sum((x - mean(x))**2) / real(n - 1, dp)
        else
            v = 0.0_dp
        end if
    end function variance

    pure function sd(x) result(s)
        ! Sample standard deviation of x (denominator n-1).
        real(dp), intent(in) :: x(:)  ! input array
        real(dp) :: s
        s = sqrt(variance(x))
    end function sd

    pure real(dp) function covariance(x, y) result(cov)
        ! Sample covariance of x and y.
        real(dp), intent(in) :: x(:), y(:)
        integer :: n

        if (size(x) /= size(y)) error stop "covariance: array sizes differ"
        n = size(x)
        if (n > 1) then
            cov = sum((x - mean(x))*(y - mean(y))) / real(n - 1, dp)
        else
            cov = 0.0_dp
        end if
    end function covariance

    pure real(dp) function correlation(x, y) result(rho)
        ! Pearson correlation of x and y.
        real(dp), intent(in) :: x(:), y(:)
        real(dp) :: xbar, ybar, xden, yden, denom

        if (size(x) /= size(y)) error stop "correlation: array sizes differ"
        if (size(x) < 2) then
            rho = 0.0_dp
            return
        end if
        xbar = mean(x)
        ybar = mean(y)
        xden = sum((x - xbar)**2)
        yden = sum((y - ybar)**2)
        denom = sqrt(xden*yden)
        if (denom <= 0.0_dp) then
            rho = 0.0_dp
        else
            rho = sum((x - xbar)*(y - ybar)) / denom
        end if
    end function correlation

    pure subroutine covariance_matrix(x, cov)
        ! Sample covariance matrix of columns of x.
        real(dp), intent(in) :: x(:, :)
        real(dp), intent(out) :: cov(:, :)
        integer :: i, j, nvar

        nvar = size(x, 2)
        if (size(cov, 1) /= nvar .or. size(cov, 2) /= nvar) then
            error stop "covariance_matrix: output shape differs"
        end if
        do i = 1, nvar
            do j = i, nvar
                cov(i, j) = covariance(x(:, i), x(:, j))
                cov(j, i) = cov(i, j)
            end do
        end do
    end subroutine covariance_matrix

    pure subroutine correlation_matrix(x, corr)
        ! Pearson correlation matrix of columns of x.
        real(dp), intent(in) :: x(:, :)
        real(dp), intent(out) :: corr(:, :)
        integer :: i, j, nvar

        nvar = size(x, 2)
        if (size(corr, 1) /= nvar .or. size(corr, 2) /= nvar) then
            error stop "correlation_matrix: output shape differs"
        end if
        do i = 1, nvar
            do j = i, nvar
                corr(i, j) = correlation(x(:, i), x(:, j))
                corr(j, i) = corr(i, j)
            end do
        end do
    end subroutine correlation_matrix

    subroutine demean_first(x, n)
        ! Subtract mean(x(1:n)) from all elements of x.
        real(dp), intent(inout) :: x(:)
        integer, intent(in) :: n

        if (n < 1 .or. n > size(x)) error stop "demean_first: invalid n"
        x = x - mean(x(1:n))
    end subroutine demean_first

    subroutine nonnegative_least_squares(x, y, beta, sse, max_iter, tol, niter)
        ! Fit y ~= x*beta by coordinate descent subject to beta >= 0.
        real(dp), intent(in) :: x(:, :), y(:)
        real(dp), intent(out) :: beta(:), sse
        integer, intent(in), optional :: max_iter
        real(dp), intent(in), optional :: tol
        integer, intent(out), optional :: niter
        real(dp), allocatable :: resid(:), xnorm2(:)
        real(dp) :: old_beta, new_beta, max_change, eps
        integer :: iter, j, nobs, nvar, iter_max, actual_iter

        nobs = size(x, 1)
        nvar = size(x, 2)
        if (size(y) /= nobs) error stop "nonnegative_least_squares: y size differs"
        if (size(beta) /= nvar) error stop "nonnegative_least_squares: beta size differs"
        iter_max = 1000
        if (present(max_iter)) iter_max = max_iter
        eps = 1.0e-10_dp
        if (present(tol)) eps = tol
        if (iter_max < 1) error stop "nonnegative_least_squares: max_iter must be positive"

        allocate(resid(nobs), xnorm2(nvar))
        beta = 0.0_dp
        resid = y
        do j = 1, nvar
            xnorm2(j) = sum(x(:, j)**2)
        end do

        actual_iter = 0
        do iter = 1, iter_max
            actual_iter = iter
            max_change = 0.0_dp
            do j = 1, nvar
                if (xnorm2(j) <= 0.0_dp) cycle
                old_beta = beta(j)
                resid = resid + x(:, j)*old_beta
                new_beta = max(0.0_dp, sum(x(:, j)*resid) / xnorm2(j))
                beta(j) = new_beta
                resid = resid - x(:, j)*new_beta
                max_change = max(max_change, abs(new_beta - old_beta))
            end do
            if (max_change <= eps*max(1.0_dp, maxval(beta))) exit
        end do

        sse = sum(resid**2)
        if (present(niter)) niter = actual_iter
        deallocate(resid, xnorm2)
    end subroutine nonnegative_least_squares

    function median(x) result(med)
        ! Sample median of x.
        real(dp), intent(in) :: x(:)
        real(dp) :: med
        real(dp), allocatable :: sorted(:)
        integer :: n

        n = size(x)
        if (n < 1) then
            med = 0.0_dp
            return
        end if
        allocate(sorted(n))
        sorted = x
        call sort_real(sorted)
        if (mod(n, 2) == 1) then
            med = sorted((n + 1) / 2)
        else
            med = 0.5_dp * (sorted(n / 2) + sorted(n / 2 + 1))
        end if
    end function median

    pure function skewness(x) result(skew)
        ! Moment skewness using central moments with denominator n.
        real(dp), intent(in) :: x(:)
        real(dp) :: skew, mu, m2, m3
        integer :: n

        n = size(x)
        if (n < 1) then
            skew = 0.0_dp
            return
        end if
        mu = mean(x)
        m2 = sum((x - mu)**2) / real(n, dp)
        if (m2 <= 0.0_dp) then
            skew = 0.0_dp
        else
            m3 = sum((x - mu)**3) / real(n, dp)
            skew = m3 / m2**1.5_dp
        end if
    end function skewness

    pure function excess_kurtosis(x) result(ekurt)
        ! Moment excess kurtosis using central moments with denominator n.
        real(dp), intent(in) :: x(:)
        real(dp) :: ekurt, mu, m2, m4
        integer :: n

        n = size(x)
        if (n < 1) then
            ekurt = 0.0_dp
            return
        end if
        mu = mean(x)
        m2 = sum((x - mu)**2) / real(n, dp)
        if (m2 <= 0.0_dp) then
            ekurt = 0.0_dp
        else
            m4 = sum((x - mu)**4) / real(n, dp)
            ekurt = m4 / m2**2 - 3.0_dp
        end if
    end function excess_kurtosis

    pure subroutine sort_real(x)
        ! Sort x in ascending order in place.
        real(dp), intent(inout) :: x(:)
        real(dp) :: key
        integer :: i, j

        do i = 2, size(x)
            key = x(i)
            j = i - 1
            do while (j >= 1)
                if (x(j) <= key) exit
                x(j+1) = x(j)
                j = j - 1
            end do
            x(j+1) = key
        end do
    end subroutine sort_real

    pure subroutine autocorr(x, acf)
        ! Autocorrelations of x at lags 1:size(acf).
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: acf(:)
        real(dp) :: mu, denom, num
        integer :: lag, i, n

        n = size(x)
        if (n < 1) then
            acf = 0.0_dp
            return
        end if
        mu = sum(x) / real(n, dp)
        denom = sum((x - mu)**2)
        if (denom <= 0.0_dp) then
            acf = 0.0_dp
            return
        end if
        do lag = 1, size(acf)
            num = 0.0_dp
            do i = lag + 1, n
                num = num + (x(i) - mu)*(x(i - lag) - mu)
            end do
            acf(lag) = num / denom
        end do
    end subroutine autocorr

    subroutine print_acf_table(title, row_label, count_label, row_values, counts, acf)
        ! Print a table of autocorrelations with one row per row_values element.
        character(len=*), intent(in) :: title, row_label, count_label
        integer, intent(in) :: row_values(:), counts(:)
        real(dp), intent(in) :: acf(:, :)
        integer :: i, lag, nlag
        character(len=16) :: lag_label
        character(len=:), allocatable :: sep

        if (size(row_values) /= size(counts)) error stop "print_acf_table: row_values/counts size mismatch"
        if (size(row_values) /= size(acf, 1)) error stop "print_acf_table: row_values/acf row size mismatch"

        nlag = size(acf, 2)
        sep = repeat("-", max(28 + 9*nlag, 44))

        print '(A)', trim(title)
        print '(A)', sep
        write(*,'(A8,1X,A10)', advance='no') row_label, count_label
        do lag = 1, nlag
            write(lag_label, '(A,I0)') "acf_", lag
            write(*,'(1X,A8)', advance='no') trim(lag_label)
        end do
        write(*,*)
        print '(A)', sep
        do i = 1, size(row_values)
            write(*,'(I8,1X,I10)', advance='no') row_values(i), counts(i)
            do lag = 1, nlag
                write(*,'(1X,F8.4)', advance='no') acf(i, lag)
            end do
            write(*,*)
        end do
        print '(A)', sep
        print '(A)', ""
    end subroutine print_acf_table

    subroutine column_summary_stats(values, med, avg, sdev, skew, ekurt, xmin, xmax)
        ! Compute median, mean, sd, skewness, excess kurtosis, min, and max for each column.
        real(dp), intent(in) :: values(:, :)
        real(dp), allocatable, intent(out) :: med(:), avg(:), sdev(:), skew(:), ekurt(:), xmin(:), xmax(:)
        integer :: j, ncols

        ncols = size(values, 2)
        allocate(med(ncols), avg(ncols), sdev(ncols), skew(ncols), ekurt(ncols), xmin(ncols), xmax(ncols))
        do j = 1, ncols
            med(j) = median(values(:, j))
            avg(j) = mean(values(:, j))
            sdev(j) = sd(values(:, j))
            skew(j) = skewness(values(:, j))
            ekurt(j) = excess_kurtosis(values(:, j))
            xmin(j) = minval(values(:, j))
            xmax(j) = maxval(values(:, j))
        end do
    end subroutine column_summary_stats

    subroutine column_robust_tail_ratios(values, tail_ratio)
        ! Compute (q95 - q05)/(q75 - q25) for each column.
        real(dp), intent(in) :: values(:, :)
        real(dp), allocatable, intent(out) :: tail_ratio(:)
        real(dp), allocatable :: sorted(:)
        real(dp) :: q05, q25, q75, q95, iqr
        integer :: j, nobs, ncols

        nobs = size(values, 1)
        ncols = size(values, 2)
        allocate(tail_ratio(ncols), sorted(nobs))
        do j = 1, ncols
            sorted = values(:, j)
            call sort_real(sorted)
            q05 = quantile_sorted(sorted, 0.05_dp)
            q25 = quantile_sorted(sorted, 0.25_dp)
            q75 = quantile_sorted(sorted, 0.75_dp)
            q95 = quantile_sorted(sorted, 0.95_dp)
            iqr = q75 - q25
            if (iqr > 0.0_dp) then
                tail_ratio(j) = (q95 - q05) / iqr
            else
                tail_ratio(j) = 2.438_dp
            end if
        end do
    end subroutine column_robust_tail_ratios

    pure function quantile_sorted(x, p) result(q)
        ! Linear-interpolated sample quantile of sorted x for 0 <= p <= 1.
        real(dp), intent(in) :: x(:), p
        real(dp) :: q, h, frac
        integer :: n, lo, hi

        n = size(x)
        if (n < 1) then
            q = 0.0_dp
            return
        end if
        h = 1.0_dp + min(max(p, 0.0_dp), 1.0_dp) * real(n - 1, dp)
        lo = int(floor(h))
        hi = min(lo + 1, n)
        frac = h - real(lo, dp)
        q = (1.0_dp - frac) * x(lo) + frac * x(hi)
    end function quantile_sorted

    subroutine print_column_summary_stats(title, col_names, med, avg, sdev, skew, ekurt, xmin, xmax)
        ! Print precomputed median, mean, sd, skewness, excess kurtosis, min, and max for each column.
        character(len=*), intent(in) :: title
        character(len=*), intent(in) :: col_names(:)
        real(dp), intent(in) :: med(:), avg(:), sdev(:), skew(:), ekurt(:), xmin(:), xmax(:)
        integer :: j, ncols

        ncols = size(col_names)
        if (size(med) /= ncols .or. size(avg) /= ncols .or. size(sdev) /= ncols .or. &
            size(skew) /= ncols .or. size(ekurt) /= ncols .or. size(xmin) /= ncols .or. &
            size(xmax) /= ncols) error stop "print_column_summary_stats: summary size mismatch"

        print '(A)', trim(title)
        print '(A)', repeat("-", 106)
        print '(A12,1X,A12,1X,A12,1X,A12,1X,A12,1X,A12,1X,A12,1X,A12)', &
            "Series", "median", "mean", "sd", "skew", "ex_kurt", "min", "max"
        print '(A)', repeat("-", 106)
        do j = 1, ncols
            print '(A12,1X,ES12.4,1X,ES12.4,1X,ES12.4,1X,F12.4,1X,F12.4,1X,ES12.4,1X,ES12.4)', &
                trim(col_names(j)), med(j), avg(j), sdev(j), skew(j), ekurt(j), xmin(j), xmax(j)
        end do
        print '(A,/)', repeat("-", 106)
    end subroutine print_column_summary_stats

end module stats_mod

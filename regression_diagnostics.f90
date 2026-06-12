! Reusable diagnostics for regression-style models.
!
! active_nnls_tstats: approximate active-set t-statistics for NNLS coefficients
! predictor_correlations: correlations of matrix columns with a response
! rank_nonzero_by_tstat: nonzero coefficient order by descending t-statistic
! print_response_predictor_corr_matrix: labeled response/predictor correlation matrix

module regression_diagnostics_mod
    use kind_mod, only: dp
    use linalg_mod, only: chol_factor, chol_solve_vec
    use stats_mod, only: correlation
    implicit none
    private
    public :: active_nnls_tstats, predictor_correlations, rank_nonzero_by_tstat
    public :: print_response_predictor_corr_matrix

contains

    subroutine active_nnls_tstats(xreg, y, beta, sse, tstat)
        ! Compute approximate OLS-style t-statistics on the active NNLS set.
        real(dp), intent(in) :: xreg(:, :), y(:), beta(:), sse
        real(dp), intent(out) :: tstat(:)
        real(dp), allocatable :: xa(:, :), xtx(:, :), L(:, :), rhs(:), sol(:)
        real(dp) :: sigma2, ridge
        integer, allocatable :: active(:)
        integer :: i, j, k, nactive, nobs
        logical :: ok

        if (size(tstat) /= size(beta)) error stop "active_nnls_tstats: tstat size differs"
        tstat = 0.0_dp
        nobs = size(y)
        nactive = count(beta > 1.0e-8_dp)
        if (nactive < 1 .or. nobs <= nactive) return
        allocate(active(nactive))
        k = 0
        do i = 1, size(beta)
            if (beta(i) > 1.0e-8_dp) then
                k = k + 1
                active(k) = i
            end if
        end do

        allocate(xa(nobs, nactive), xtx(nactive, nactive), L(nactive, nactive), rhs(nactive), sol(nactive))
        do j = 1, nactive
            xa(:, j) = xreg(:, active(j))
        end do
        xtx = matmul(transpose(xa), xa)
        call chol_factor(xtx, nactive, L, ok)
        if (.not. ok) then
            ridge = max(1.0e-10_dp*maxval(abs(xtx)), 1.0e-12_dp)
            do i = 1, nactive
                xtx(i, i) = xtx(i, i) + ridge
            end do
            call chol_factor(xtx, nactive, L, ok)
            if (.not. ok) then
                deallocate(active, xa, xtx, L, rhs, sol)
                return
            end if
        end if

        sigma2 = sse / real(nobs - nactive, dp)
        do j = 1, nactive
            rhs = 0.0_dp
            rhs(j) = 1.0_dp
            call chol_solve_vec(L, nactive, rhs, sol)
            if (sol(j) > 0.0_dp .and. sigma2 > 0.0_dp) then
                tstat(active(j)) = beta(active(j)) / sqrt(sigma2*sol(j))
            end if
        end do
        deallocate(active, xa, xtx, L, rhs, sol)
    end subroutine active_nnls_tstats

    subroutine predictor_correlations(xreg, y, pred_corr)
        ! Compute the correlation of each xreg column with response y.
        real(dp), intent(in) :: xreg(:, :), y(:)
        real(dp), intent(out) :: pred_corr(:)
        integer :: j

        if (size(pred_corr) /= size(xreg, 2)) error stop "predictor_correlations: size mismatch"
        do j = 1, size(xreg, 2)
            pred_corr(j) = correlation(xreg(:, j), y)
        end do
    end subroutine predictor_correlations

    subroutine rank_nonzero_by_tstat(beta, tstat, order)
        ! Return nonzero beta indices sorted by descending t-statistic; unused entries are zero.
        real(dp), intent(in) :: beta(:), tstat(:)
        integer, intent(out) :: order(:)
        logical, allocatable :: used(:)
        real(dp) :: best_t
        integer :: i, k, best

        if (size(beta) /= size(tstat) .or. size(order) /= size(beta)) then
            error stop "rank_nonzero_by_tstat: size mismatch"
        end if
        allocate(used(size(beta)))
        used = .false.
        order = 0
        do k = 1, size(beta)
            best = 0
            best_t = -huge(1.0_dp)
            do i = 1, size(beta)
                if (used(i) .or. beta(i) <= 1.0e-8_dp) cycle
                if (tstat(i) > best_t) then
                    best_t = tstat(i)
                    best = i
                end if
            end do
            if (best == 0) exit
            order(k) = best
            used(best) = .true.
        end do
        deallocate(used)
    end subroutine rank_nonzero_by_tstat

    subroutine print_response_predictor_corr_matrix(title, response_label, predictor_labels, y, x)
        ! Print a compact correlation matrix for one response and several predictors.
        character(len=*), intent(in) :: title, response_label, predictor_labels(:)
        real(dp), intent(in) :: y(:), x(:, :)
        character(len=32), allocatable :: labels(:)
        character(len=6), allocatable :: codes(:)
        real(dp), allocatable :: z(:, :)
        integer :: i, j, ncol

        if (size(predictor_labels) /= size(x, 2)) then
            error stop "print_response_predictor_corr_matrix: label count differs"
        end if
        if (size(y) /= size(x, 1)) error stop "print_response_predictor_corr_matrix: row count differs"
        if (size(x, 2) < 1) return

        ncol = size(x, 2) + 1
        allocate(labels(ncol), codes(ncol), z(size(y), ncol))
        labels(1) = adjustl(response_label)
        codes(1) = "Y"
        z(:, 1) = y
        do j = 1, size(x, 2)
            write (codes(j + 1), '("P",I2.2)') j
            labels(j + 1) = adjustl(predictor_labels(j))
            z(:, j + 1) = x(:, j)
        end do

        print '(A)', ""
        print '(A)', trim(title)
        print '(A)', "Retained predictor labels"
        print '(A)', "-------------------------------------------"
        do i = 1, ncol
            print '(A6,1X,A32)', codes(i), labels(i)
        end do
        print '(A)', "-------------------------------------------"
        print '(A)', "------" // repeat("---------", ncol)
        write (*, '(A6)', advance='no') ""
        do j = 1, ncol
            write (*, '(1X,A8)', advance='no') codes(j)
        end do
        print '(A)', ""
        print '(A)', "------" // repeat("---------", ncol)
        do i = 1, ncol
            write (*, '(A6)', advance='no') codes(i)
            do j = 1, ncol
                write (*, '(1X,F8.4)', advance='no') correlation(z(:, i), z(:, j))
            end do
            print '(A)', ""
        end do
        print '(A)', "------" // repeat("---------", ncol)
        deallocate(labels, codes, z)
    end subroutine print_response_predictor_corr_matrix

end module regression_diagnostics_mod

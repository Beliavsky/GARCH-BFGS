! Reporting helper for intraday cross-asset covariance/correlation tables.

module intraday_correlation_report_mod
    use kind_mod, only: dp
    use market_data_mod, only: ohlcv_series_t, default_session_start_seconds
    use intraday_returns_mod, only: return_series_t, returns_at_frequency, aligned_return_matrix, &
                                    realized_vol_ann_pct, frequency_label
    use matrix_print_mod, only: print_real_matrix, print_vol_table
    use stats_mod, only: covariance_matrix, correlation_matrix
    implicit none
    private

    public :: print_frequency_correlation_report
    public :: print_correlation_anomaly_report

contains

    ! Print realized volatilities and covariance/correlation matrices for one frequency.
    subroutine print_frequency_correlation_report(target_seconds, source_seconds, labels, bars, regular_session, &
                                                  corr_out, nobs_out)
        integer, intent(in) :: target_seconds, source_seconds
        character(len=*), intent(in) :: labels(:)
        type(ohlcv_series_t), intent(in) :: bars(:)
        logical, intent(in) :: regular_session
        real(dp), intent(out), optional :: corr_out(:, :)
        integer, intent(out), optional :: nobs_out
        type(return_series_t), allocatable :: returns(:)
        real(dp), allocatable :: x(:, :), cov(:, :), corr(:, :), vol_ann(:), rv(:)
        integer :: nasset, nobs, i

        nasset = size(labels)
        if (present(corr_out)) corr_out = 0.0_dp
        if (present(nobs_out)) nobs_out = 0
        allocate(returns(nasset))
        do i = 1, nasset
            if (regular_session) then
                call returns_at_frequency(bars(i), source_seconds, target_seconds, returns(i), &
                                          session_start_seconds=default_session_start_seconds)
            else
                call returns_at_frequency(bars(i), source_seconds, target_seconds, returns(i), session_start_seconds=0)
            end if
        end do
        call aligned_return_matrix(returns, x)
        nobs = size(x, 1)
        if (present(nobs_out)) nobs_out = nobs
        if (nobs < 2) then
            print '(A,I0,A)', "Skipping ", target_seconds, " seconds: fewer than 2 aligned returns."
            return
        end if
        allocate(cov(nasset, nasset), corr(nasset, nasset), vol_ann(nasset), rv(nasset))
        call covariance_matrix(x, cov)
        call correlation_matrix(x, corr)
        if (present(corr_out)) corr_out = corr
        do i = 1, nasset
            rv(i) = sum(x(:, i)**2)
            vol_ann(i) = realized_vol_ann_pct(x(:, i), target_seconds, regular_session=regular_session)
        end do

        print '(A)', repeat("-", 96)
        print '(A,I0,A,A,A,I0)', "Frequency: ", target_seconds, " seconds (", frequency_label(target_seconds), &
              ")  aligned returns: ", nobs
        call print_vol_table(labels, rv, vol_ann)
        call print_real_matrix("Covariance matrix", labels, cov)
        call print_real_matrix("Correlation matrix", labels, corr, value_format="(1X,F6.3)", value_width=6)
        deallocate(returns, x, cov, corr, vol_ann, rv)
    end subroutine print_frequency_correlation_report

    ! Compare lower-frequency correlations to the first/high-frequency matrix.
    subroutine print_correlation_anomaly_report(labels, freq_labels, corr_by_freq, nobs_by_freq, &
                                                z_threshold, max_pairs)
        character(len=*), intent(in) :: labels(:), freq_labels(:)
        real(dp), intent(in) :: corr_by_freq(:, :, :)
        integer, intent(in) :: nobs_by_freq(:)
        real(dp), intent(in), optional :: z_threshold
        integer, intent(in), optional :: max_pairs
        integer :: nf, na, f, i, j, npair, nflag, ntop, max_top
        real(dp) :: threshold, expected, p_two, z, diff, se, rho0, rho1
        real(dp), allocatable :: top_abs_z(:), top_z(:), top_diff(:), top_rho0(:), top_rho1(:)
        integer, allocatable :: top_i(:), top_j(:)

        na = size(labels)
        nf = size(freq_labels)
        if (size(corr_by_freq, 1) /= na .or. size(corr_by_freq, 2) /= na .or. size(corr_by_freq, 3) /= nf) then
            error stop "print_correlation_anomaly_report: correlation array shape differs"
        end if
        if (size(nobs_by_freq) /= nf) error stop "print_correlation_anomaly_report: nobs size differs"
        if (nf < 2 .or. na < 2) return
        if (nobs_by_freq(1) < 4) return

        threshold = 3.0_dp
        if (present(z_threshold)) threshold = z_threshold
        max_top = 20
        if (present(max_pairs)) max_top = max_pairs
        max_top = max(1, max_top)
        npair = na*(na - 1)/2
        p_two = 2.0_dp*(1.0_dp - normal_cdf(threshold))
        expected = real(npair, dp)*p_two

        allocate(top_abs_z(max_top), top_z(max_top), top_diff(max_top), top_rho0(max_top), top_rho1(max_top))
        allocate(top_i(max_top), top_j(max_top))

        print '(A)', ""
        print '(A)', "Correlation changes versus highest frequency"
        print '(A,F6.2,A)', "Fisher z threshold: ", threshold, "  (IID expected counts use pairwise tests)"
        print '(A)', "--------------------------------------------------------------------------------"
        print '(A16,1X,A16,1X,A8,1X,A8,1X,A12,1X,A12)', "base_freq", "compare_freq", "n_base", "n_comp", &
              "flags", "expected"
        print '(A)', "--------------------------------------------------------------------------------"
        do f = 2, nf
            if (nobs_by_freq(f) < 4) cycle
            nflag = 0
            ntop = 0
            top_abs_z = -1.0_dp
            do i = 1, na - 1
                do j = i + 1, na
                    rho0 = clamp_corr(corr_by_freq(i, j, 1))
                    rho1 = clamp_corr(corr_by_freq(i, j, f))
                    se = sqrt(1.0_dp/real(nobs_by_freq(1) - 3, dp) + 1.0_dp/real(nobs_by_freq(f) - 3, dp))
                    z = (fisher_z(rho1) - fisher_z(rho0)) / se
                    diff = rho1 - rho0
                    if (abs(z) >= threshold) nflag = nflag + 1
                    call update_top_pairs(abs(z), z, diff, rho0, rho1, i, j, ntop, top_abs_z, top_z, &
                                          top_diff, top_rho0, top_rho1, top_i, top_j)
                end do
            end do
            print '(A16,1X,A16,1X,I8,1X,I8,1X,I12,1X,F12.3)', freq_labels(1), freq_labels(f), &
                  nobs_by_freq(1), nobs_by_freq(f), nflag, expected
            call print_top_pairs(labels, freq_labels(1), freq_labels(f), ntop, top_z, top_diff, &
                                 top_rho0, top_rho1, top_i, top_j)
        end do
        print '(A)', "--------------------------------------------------------------------------------"
        deallocate(top_abs_z, top_z, top_diff, top_rho0, top_rho1, top_i, top_j)
    end subroutine print_correlation_anomaly_report

    ! Insert one pair into the top-|z| list.
    subroutine update_top_pairs(abs_z, z, diff, rho0, rho1, i, j, ntop, top_abs_z, top_z, top_diff, &
                                top_rho0, top_rho1, top_i, top_j)
        real(dp), intent(in) :: abs_z, z, diff, rho0, rho1
        integer, intent(in) :: i, j
        integer, intent(inout) :: ntop
        real(dp), intent(inout) :: top_abs_z(:), top_z(:), top_diff(:), top_rho0(:), top_rho1(:)
        integer, intent(inout) :: top_i(:), top_j(:)
        integer :: k, pos, nmax

        nmax = size(top_abs_z)
        if (ntop < nmax) then
            ntop = ntop + 1
            pos = ntop
        else
            pos = minloc(top_abs_z, dim=1)
            if (abs_z <= top_abs_z(pos)) return
        end if
        top_abs_z(pos) = abs_z
        top_z(pos) = z
        top_diff(pos) = diff
        top_rho0(pos) = rho0
        top_rho1(pos) = rho1
        top_i(pos) = i
        top_j(pos) = j
        do k = pos, 2, -1
            if (top_abs_z(k) <= top_abs_z(k - 1)) exit
            call swap_top(k, k - 1, top_abs_z, top_z, top_diff, top_rho0, top_rho1, top_i, top_j)
        end do
    end subroutine update_top_pairs

    ! Swap two rows of the top-pair buffers.
    subroutine swap_top(a, b, top_abs_z, top_z, top_diff, top_rho0, top_rho1, top_i, top_j)
        integer, intent(in) :: a, b
        real(dp), intent(inout) :: top_abs_z(:), top_z(:), top_diff(:), top_rho0(:), top_rho1(:)
        integer, intent(inout) :: top_i(:), top_j(:)
        real(dp) :: xr
        integer :: xi

        xr = top_abs_z(a); top_abs_z(a) = top_abs_z(b); top_abs_z(b) = xr
        xr = top_z(a); top_z(a) = top_z(b); top_z(b) = xr
        xr = top_diff(a); top_diff(a) = top_diff(b); top_diff(b) = xr
        xr = top_rho0(a); top_rho0(a) = top_rho0(b); top_rho0(b) = xr
        xr = top_rho1(a); top_rho1(a) = top_rho1(b); top_rho1(b) = xr
        xi = top_i(a); top_i(a) = top_i(b); top_i(b) = xi
        xi = top_j(a); top_j(a) = top_j(b); top_j(b) = xi
    end subroutine swap_top

    ! Print top pairwise correlation changes for one frequency comparison.
    subroutine print_top_pairs(labels, base_label, comp_label, ntop, top_z, top_diff, top_rho0, top_rho1, top_i, top_j)
        character(len=*), intent(in) :: labels(:), base_label, comp_label
        integer, intent(in) :: ntop, top_i(:), top_j(:)
        real(dp), intent(in) :: top_z(:), top_diff(:), top_rho0(:), top_rho1(:)
        integer :: k

        print '(A,A,A,A)', "Top pair changes: ", trim(base_label), " vs ", trim(comp_label)
        print '(A16,1X,A16,1X,A8,1X,A8,1X,A8,1X,A8)', "asset_1", "asset_2", "rho_base", "rho_comp", &
              "diff", "z_diff"
        do k = 1, ntop
            print '(A16,1X,A16,1X,F8.3,1X,F8.3,1X,F8.3,1X,F8.2)', labels(top_i(k)), labels(top_j(k)), &
                  top_rho0(k), top_rho1(k), top_diff(k), top_z(k)
        end do
    end subroutine print_top_pairs

    pure real(dp) function clamp_corr(rho) result(x)
        real(dp), intent(in) :: rho

        x = min(max(rho, -0.999999_dp), 0.999999_dp)
    end function clamp_corr

    pure real(dp) function fisher_z(rho) result(z)
        real(dp), intent(in) :: rho

        z = 0.5_dp*log((1.0_dp + rho) / (1.0_dp - rho))
    end function fisher_z

    pure real(dp) function normal_cdf(x) result(p)
        real(dp), intent(in) :: x

        p = 0.5_dp*(1.0_dp + erf(x / sqrt(2.0_dp)))
    end function normal_cdf

end module intraday_correlation_report_mod

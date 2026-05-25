module garch_forecast_mod
    use kind_mod,        only: dp
    use csv_mod,         only: date_label
    use stats_mod,       only: sort_real
    use garch_types_mod, only: garch_params_t, garch_fit_result_t
    use garch_fit_mod,   only: riskmetrics2006_variance, ewma_variance, figarch_variance, fi_nagarch_variance, &
                               symm_garch_pq_variance, nagarch_pq_variance, &
                               csgarch_variance, tgarch_variance, avgarch_variance, &
                               garch_skew_kurt, garch_pq_skew_kurt, qgarch_skew_kurt, &
                               figarch_skew_kurt, fi_nagarch_skew_kurt, &
                               nagarch_skew_kurt, nagarch_pq_skew_kurt, csgarch_skew_kurt, gjr_skew_kurt, &
                               aparch_skew_kurt, harch_skew_kurt, tgarch_skew_kurt, avgarch_skew_kurt, riskmetrics2006_skew_kurt, &
                               ewma_skew_kurt, &
                               midas_hyperbolic_skew_kurt, midas_hyperbolic_asym_skew_kurt, &
                               symm_garch_persist, symm_garch_pq_persist, qgarch_persist, &
                               figarch_persist, fi_nagarch_persist, &
                               nagarch_persist, nagarch_pq_persist, gjr_persist, fgarch_twist_persist, &
                               csgarch_persist, aparch_persist, harch_persist, tgarch_persist, avgarch_persist, riskmetrics2006_persist, &
                               ewma_persist, &
                               midas_hyperbolic_persist, midas_hyperbolic_asym_persist, &
                               aparch_mean_variance, qgarch_mean_variance
    implicit none
    private

    public :: model_vol_forecast, finalize_return_garch_fit
    public :: vol_forecast_stats_t, summarize_vol_forecast, print_vol_forecast_table

    type :: vol_forecast_stats_t
        real(dp) :: median = 0.0_dp
        real(dp) :: mean = 0.0_dp
        real(dp) :: sd = 0.0_dp
        real(dp) :: ekurt = 0.0_dp
        real(dp) :: vmin = 0.0_dp
        real(dp) :: vmax = 0.0_dp
        real(dp) :: first = 0.0_dp
        real(dp) :: last = 0.0_dp
        integer :: date_min = 0
        integer :: date_max = 0
    end type vol_forecast_stats_t

contains

    subroutine finalize_return_garch_fit(model_name, y, params, fopt, niter, converged, trading_days, result, &
                                         vol_ann_override, skew_override, ekurt_override)
        character(len=*), intent(in) :: model_name
        real(dp), intent(in) :: y(:), fopt, trading_days
        type(garch_params_t), intent(in) :: params
        integer, intent(in) :: niter
        logical, intent(in) :: converged
        type(garch_fit_result_t), intent(out) :: result
        real(dp), intent(in), optional :: vol_ann_override, skew_override, ekurt_override
        real(dp) :: h_unc
        real(dp), allocatable :: variance(:)
        integer :: nobs

        nobs = size(y)
        result = garch_fit_result_t()
        result%model = trim(model_name)
        result%params = params
        result%niter = niter
        result%converged = converged
        result%logl = -real(nobs, dp) * fopt

        select case (trim(model_name))
        case ("EWMA")
            result%persist = ewma_persist(params)
            allocate(variance(nobs))
            call ewma_variance(y, params, variance)
            h_unc = sum(variance) / real(nobs, dp)
            deallocate(variance)
            result%nparam = 1
            call ewma_skew_kurt(y, params, result%skew, result%ekurt)
        case ("SYMM_GARCH")
            result%persist = symm_garch_persist(params)
            h_unc = params%omega / max(1.0_dp - result%persist, 1.0e-8_dp)
            result%nparam = 3
            call garch_skew_kurt(y, params, result%skew, result%ekurt)
        case ("SYMM_GARCH_2_1", "SYMM_GARCH_1_2", "SYMM_GARCH_2_2")
            result%persist = symm_garch_pq_persist(params)
            allocate(variance(nobs))
            call symm_garch_pq_variance(y, params, variance)
            h_unc = sum(variance) / real(nobs, dp)
            deallocate(variance)
            result%nparam = 1 + size(params%alpha_lags) + size(params%beta_lags)
            call garch_pq_skew_kurt(y, params, result%skew, result%ekurt)
        case ("QGARCH")
            result%persist = qgarch_persist(params)
            h_unc = qgarch_mean_variance(params)
            result%nparam = 4
            call qgarch_skew_kurt(y, params, result%skew, result%ekurt)
        case ("CSGARCH")
            result%persist = csgarch_persist(params)
            allocate(variance(nobs))
            call csgarch_variance(y, params, variance)
            h_unc = sum(variance) / real(nobs, dp)
            deallocate(variance)
            result%nparam = 5
            call csgarch_skew_kurt(y, params, result%skew, result%ekurt)
        case ("FIGARCH")
            result%persist = figarch_persist(params)
            allocate(variance(nobs))
            call figarch_variance(y, params, variance)
            h_unc = sum(variance) / real(nobs, dp)
            deallocate(variance)
            result%nparam = 4
            call figarch_skew_kurt(y, params, result%skew, result%ekurt)
        case ("FI_NAGARCH")
            result%persist = fi_nagarch_persist(params)
            allocate(variance(nobs))
            call fi_nagarch_variance(y, params, variance)
            h_unc = sum(variance) / real(nobs, dp)
            deallocate(variance)
            result%nparam = 5
            call fi_nagarch_skew_kurt(y, params, result%skew, result%ekurt)
        case ("NAGARCH")
            result%persist = nagarch_persist(params)
            h_unc = params%omega / max(1.0_dp - result%persist, 1.0e-8_dp)
            result%nparam = 4
            call nagarch_skew_kurt(y, params, result%skew, result%ekurt)
        case ("NAGARCH_2_1", "NAGARCH_1_2", "NAGARCH_2_2")
            result%persist = nagarch_pq_persist(params)
            allocate(variance(nobs))
            call nagarch_pq_variance(y, params, variance)
            h_unc = sum(variance) / real(nobs, dp)
            deallocate(variance)
            result%nparam = 2 + size(params%alpha_lags) + size(params%beta_lags)
            call nagarch_pq_skew_kurt(y, params, result%skew, result%ekurt)
        case ("GJR_GARCH")
            result%persist = gjr_persist(params)
            h_unc = params%omega / max(1.0_dp - result%persist, 1.0e-8_dp)
            result%nparam = 4
            call gjr_skew_kurt(y, params, result%skew, result%ekurt)
        case ("FGTWIST")
            result%persist = fgarch_twist_persist(params)
            h_unc = params%omega / max(1.0_dp - result%persist, 1.0e-8_dp)
            result%nparam = 5
        case ("APARCH")
            result%persist = aparch_persist(params)
            h_unc = aparch_mean_variance(y, params)
            result%nparam = 5
            call aparch_skew_kurt(y, params, result%skew, result%ekurt)
        case ("HARCH")
            result%persist = harch_persist(params)
            h_unc = params%omega / max(1.0_dp - result%persist, 1.0e-8_dp)
            result%nparam = 4
            call harch_skew_kurt(y, params, result%skew, result%ekurt)
        case ("TGARCH")
            result%persist = tgarch_persist(params)
            allocate(variance(nobs))
            call tgarch_variance(y, params, variance)
            h_unc = sum(variance) / real(nobs, dp)
            deallocate(variance)
            result%nparam = 4
            call tgarch_skew_kurt(y, params, result%skew, result%ekurt)
        case ("AVGARCH")
            result%persist = avgarch_persist(params)
            allocate(variance(nobs))
            call avgarch_variance(y, params, variance)
            h_unc = sum(variance) / real(nobs, dp)
            deallocate(variance)
            result%nparam = 5
            call avgarch_skew_kurt(y, params, result%skew, result%ekurt)
        case ("RM2006")
            result%persist = riskmetrics2006_persist()
            allocate(variance(nobs))
            call riskmetrics2006_variance(y, variance)
            h_unc = sum(variance) / real(nobs, dp)
            deallocate(variance)
            result%nparam = 0
            call riskmetrics2006_skew_kurt(y, result%skew, result%ekurt)
        case ("MIDASHYP")
            result%persist = midas_hyperbolic_persist(params)
            h_unc = params%omega / max(1.0_dp - result%persist, 1.0e-8_dp)
            result%nparam = 3
            call midas_hyperbolic_skew_kurt(y, params, result%skew, result%ekurt)
        case ("MIDASHYP_ASYM")
            result%persist = midas_hyperbolic_asym_persist(params)
            h_unc = params%omega / max(1.0_dp - result%persist, 1.0e-8_dp)
            result%nparam = 4
            call midas_hyperbolic_asym_skew_kurt(y, params, result%skew, result%ekurt)
        case default
            h_unc = params%omega / max(1.0_dp - result%persist, 1.0e-8_dp)
        end select

        result%vol_ann = sqrt(trading_days * h_unc) * 100.0_dp
        if (present(vol_ann_override)) result%vol_ann = vol_ann_override
        if (present(skew_override)) result%skew = skew_override
        if (present(ekurt_override)) result%ekurt = ekurt_override
        result%aic = 2.0_dp * real(result%nparam, dp) - 2.0_dp * result%logl
        result%bic = log(real(nobs, dp)) * real(result%nparam, dp) - 2.0_dp * result%logl
    end subroutine finalize_return_garch_fit

    subroutine model_vol_forecast(model_name, y, params, persist, trading_days, vol)
        character(len=*), intent(in) :: model_name
        real(dp), intent(in) :: y(:)
        type(garch_params_t), intent(in) :: params
        real(dp), intent(in) :: persist, trading_days
        real(dp), intent(out) :: vol(:)
        real(dp) :: h, sqrth, z, r, ind, q, sdel, backcast
        real(dp), allocatable :: variance(:)
        integer :: t

        if (trim(model_name) == "RM2006") then
            allocate(variance(size(y)))
            call riskmetrics2006_variance(y, variance)
            vol = sqrt(trading_days * variance) * 100.0_dp
            deallocate(variance)
            return
        else if (trim(model_name) == "FIGARCH") then
            allocate(variance(size(y)))
            call figarch_variance(y, params, variance)
            vol = sqrt(trading_days * variance) * 100.0_dp
            deallocate(variance)
            return
        else if (trim(model_name) == "FI_NAGARCH") then
            allocate(variance(size(y)))
            call fi_nagarch_variance(y, params, variance)
            vol = sqrt(trading_days * variance) * 100.0_dp
            deallocate(variance)
            return
        else if (trim(model_name) == "SYMM_GARCH_2_1" .or. trim(model_name) == "SYMM_GARCH_1_2" .or. &
            trim(model_name) == "SYMM_GARCH_2_2") then
            allocate(variance(size(y)))
            call symm_garch_pq_variance(y, params, variance)
            vol = sqrt(trading_days * variance) * 100.0_dp
            deallocate(variance)
            return
        else if (trim(model_name) == "CSGARCH") then
            allocate(variance(size(y)))
            call csgarch_variance(y, params, variance)
            vol = sqrt(trading_days * variance) * 100.0_dp
            deallocate(variance)
            return
        else if (trim(model_name) == "TGARCH") then
            allocate(variance(size(y)))
            call tgarch_variance(y, params, variance)
            vol = sqrt(trading_days * variance) * 100.0_dp
            deallocate(variance)
            return
        else if (trim(model_name) == "AVGARCH") then
            allocate(variance(size(y)))
            call avgarch_variance(y, params, variance)
            vol = sqrt(trading_days * variance) * 100.0_dp
            deallocate(variance)
            return
        else if (trim(model_name) == "NAGARCH_2_1" .or. trim(model_name) == "NAGARCH_1_2" .or. &
            trim(model_name) == "NAGARCH_2_2") then
            allocate(variance(size(y)))
            call nagarch_pq_variance(y, params, variance)
            vol = sqrt(trading_days * variance) * 100.0_dp
            deallocate(variance)
            return
        end if

        sdel = 0.0_dp
        backcast = max(sum(y**2) / real(size(y), dp), 1.0e-12_dp)
        if (trim(model_name) == "APARCH") then
            sdel = max(sum(abs(y)**params%theta) / real(size(y), dp), 1.0e-12_dp)
            h = sdel**(2.0_dp / params%theta)
        else if (trim(model_name) == "EWMA" .or. trim(model_name) == "AEWMA_NAG" .or. &
            trim(model_name) == "AEWMA_TWIST") then
            h = sum(y**2) / real(size(y), dp)
        else if (trim(model_name) == "QGARCH") then
            h = (params%omega + params%alpha*params%theta**2) / max(1.0_dp - persist, 1.0e-8_dp)
        else if (trim(model_name) == "HARCH") then
            h = params%omega / max(1.0_dp - persist, 1.0e-8_dp)
        else
            h = params%omega / max(1.0_dp - persist, 1.0e-8_dp)
        end if

        do t = 1, size(y)
            h = max(h, 1.0e-12_dp)
            sqrth = sqrt(h)
            vol(t) = sqrt(trading_days * h) * 100.0_dp
            select case (trim(model_name))
            case ("EWMA")
                h = params%beta*h + params%alpha*y(t)**2
            case ("AEWMA_NAG")
                z = y(t) / sqrth
                q = z - params%theta
                h = params%beta*h + params%alpha*params%scale*h*q**2
            case ("AEWMA_TWIST")
                z = y(t) / sqrth
                q = abs(z - params%theta) - params%twist*(z - params%theta)
                h = params%beta*h + params%alpha*params%scale*h*q**2
            case ("SYMM_GARCH")
                h = params%omega + params%alpha*y(t)**2 + params%beta*h
            case ("QGARCH")
                h = params%omega + params%alpha*(y(t) - params%theta)**2 + params%beta*h
            case ("NAGARCH")
                r = y(t) - params%theta*sqrth
                h = params%omega + params%alpha*r**2 + params%beta*h
            case ("GJR_GARCH")
                ind = merge(1.0_dp, 0.0_dp, y(t) < 0.0_dp)
                h = params%omega + (params%alpha + params%gamma*ind)*y(t)**2 + params%beta*h
            case ("FGTWIST")
                z = y(t) / sqrth
                q = abs(z - params%theta) - params%twist*(z - params%theta)
                h = params%omega + params%alpha*h*q**2 + params%beta*h
            case ("APARCH")
                q = max(abs(y(t)) - params%gamma*y(t), 1.0e-12_dp)
                sdel = params%omega + params%alpha*q**params%theta + params%beta*sdel
                h = max(sdel, 1.0e-12_dp)**(2.0_dp / params%theta)
            case ("HARCH")
                h = params%omega + params%alpha*harch_block_y(y, t, 1, backcast) + &
                    params%gamma*harch_block_y(y, t, 5, backcast) + &
                    params%beta*harch_block_y(y, t, 22, backcast)
            case ("MIDASHYP")
                h = params%omega + params%alpha*midas_hyperbolic_block_y(y, t, params%theta, backcast)
            case ("MIDASHYP_ASYM")
                h = params%omega + midas_hyperbolic_asym_block_y(y, t, params%alpha, params%gamma, params%theta, backcast)
            end select
        end do
    end subroutine model_vol_forecast

    subroutine summarize_vol_forecast(vol, return_dates, stats)
        real(dp), intent(in) :: vol(:)
        integer,  intent(in) :: return_dates(:)
        type(vol_forecast_stats_t), intent(out) :: stats
        real(dp), allocatable :: sorted(:)
        real(dp) :: rn, dz, m2, m4
        integer :: n, t, imin, imax

        n = size(vol)
        rn = real(n, dp)
        allocate(sorted(n))
        sorted = vol
        call sort_real(sorted)
        if (mod(n, 2) == 1) then
            stats%median = sorted((n + 1) / 2)
        else
            stats%median = 0.5_dp * (sorted(n / 2) + sorted(n / 2 + 1))
        end if
        stats%mean = sum(vol) / rn
        if (n > 1) then
            stats%sd = sqrt(sum((vol - stats%mean)**2) / real(n - 1, dp))
        else
            stats%sd = 0.0_dp
        end if
        m2 = 0.0_dp
        m4 = 0.0_dp
        do t = 1, n
            dz = vol(t) - stats%mean
            m2 = m2 + dz**2
            m4 = m4 + dz**4
        end do
        if (m2 > 0.0_dp) then
            stats%ekurt = rn*m4/(m2*m2) - 3.0_dp
        else
            stats%ekurt = 0.0_dp
        end if
        imin = minloc(vol, dim=1)
        imax = maxloc(vol, dim=1)
        stats%vmin = vol(imin)
        stats%vmax = vol(imax)
        stats%first = vol(1)
        stats%last = vol(n)
        stats%date_min = return_dates(imin)
        stats%date_max = return_dates(imax)
        deallocate(sorted)
    end subroutine summarize_vol_forecast

    subroutine print_vol_forecast_table(asset_names, model_names, stats, have)
        character(len=*), intent(in) :: asset_names(:), model_names(:)
        type(vol_forecast_stats_t), intent(in) :: stats(:,:)
        logical, intent(in) :: have(:,:)
        integer :: iasset, imodel

        write(*,'(/,A)') "Volatility forecast properties:"
        write(*,'(A)') "Model            Asset       median      mean        sd     ekurt       min       max     first      last date_of_min date_of_max"
        write(*,'(A)') repeat("-", 128)
        do iasset = 1, size(asset_names)
            do imodel = 1, size(model_names)
                if (.not. have(iasset,imodel)) cycle
                write(*,'(A16,1X,A9,8F10.3,2(1X,A10))') trim(model_names(imodel)), trim(asset_names(iasset)), &
                    stats(iasset,imodel)%median, stats(iasset,imodel)%mean, stats(iasset,imodel)%sd, &
                    stats(iasset,imodel)%ekurt, stats(iasset,imodel)%vmin, stats(iasset,imodel)%vmax, &
                    stats(iasset,imodel)%first, stats(iasset,imodel)%last, &
                    date_label(stats(iasset,imodel)%date_min), date_label(stats(iasset,imodel)%date_max)
            end do
        end do
    end subroutine print_vol_forecast_table

    real(dp) function midas_hyperbolic_block_y(y, t, theta, backcast)
        real(dp), intent(in) :: y(:), theta, backcast
        integer, intent(in) :: t
        real(dp) :: raw(22), sum_raw, ratio, lag_sq
        integer :: i, idx

        raw(1) = theta
        do i = 2, 22
            ratio = (real(i - 1, dp) + theta) / real(i, dp)
            raw(i) = raw(i - 1) * ratio
        end do
        sum_raw = sum(raw)
        midas_hyperbolic_block_y = 0.0_dp
        do i = 1, 22
            idx = t - i
            if (idx >= 1) then
                lag_sq = y(idx)**2
            else
                lag_sq = backcast
            end if
            midas_hyperbolic_block_y = midas_hyperbolic_block_y + raw(i)*lag_sq / sum_raw
        end do
    end function midas_hyperbolic_block_y

    real(dp) function midas_hyperbolic_asym_block_y(y, t, alpha, gamma, theta, backcast)
        real(dp), intent(in) :: y(:), alpha, gamma, theta, backcast
        integer, intent(in) :: t
        real(dp) :: raw(22), sum_raw, ratio, lag_sq, ind
        integer :: i, idx

        raw(1) = theta
        do i = 2, 22
            ratio = (real(i - 1, dp) + theta) / real(i, dp)
            raw(i) = raw(i - 1) * ratio
        end do
        sum_raw = sum(raw)
        midas_hyperbolic_asym_block_y = 0.0_dp
        do i = 1, 22
            idx = t - i
            if (idx >= 1) then
                lag_sq = y(idx)**2
                ind = merge(1.0_dp, 0.0_dp, y(idx) < 0.0_dp)
            else
                lag_sq = backcast
                ind = 0.5_dp
            end if
            midas_hyperbolic_asym_block_y = midas_hyperbolic_asym_block_y + &
                (alpha + gamma*ind)*raw(i)*lag_sq / sum_raw
        end do
    end function midas_hyperbolic_asym_block_y

    real(dp) function harch_block_y(y, t, lag, backcast)
        real(dp), intent(in) :: y(:), backcast
        integer, intent(in) :: t, lag
        integer :: i, idx

        harch_block_y = 0.0_dp
        do i = 1, lag
            idx = t - i
            if (idx >= 1) then
                harch_block_y = harch_block_y + y(idx)**2
            else
                harch_block_y = harch_block_y + backcast
            end if
        end do
        harch_block_y = harch_block_y / real(lag, dp)
    end function harch_block_y

end module garch_forecast_mod

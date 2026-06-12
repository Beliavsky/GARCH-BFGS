! Utilities for forecasting variance remaining until the market close.
!
! A diurnal curve is treated as a positive intraday variance shape. Its scale
! does not matter for fractions: remaining fractions are computed after
! normalizing the curve to sum to one across the session.

module remaining_vol_mod
    use kind_mod, only: dp
    implicit none
    private

    real(dp), parameter :: min_var = 1.0e-12_dp

    public :: normalize_diurnal_weights
    public :: cumulative_diurnal_fraction
    public :: remaining_fraction_from_bin
    public :: remaining_fraction_after_bin
    public :: remaining_session_fraction
    public :: remaining_day_fraction_from_bin
    public :: remaining_variance_from_curve
    public :: remaining_forecast_from_curve
    public :: remaining_variance_from_path
    public :: full_session_variance_from_path
    public :: annualized_vol_pct
    public :: trading_day_annualized_vol_pct
    public :: compose_expiry_variance

contains

    ! Normalize a positive diurnal variance curve so its entries sum to one.
    pure subroutine normalize_diurnal_weights(diurnal_var, weight)
        real(dp), intent(in) :: diurnal_var(:)
        real(dp), intent(out) :: weight(:)
        real(dp) :: total

        if (size(weight) /= size(diurnal_var)) then
            error stop "normalize_diurnal_weights: array sizes differ"
        end if
        if (size(diurnal_var) < 1) error stop "normalize_diurnal_weights: empty curve"
        weight = max(diurnal_var, 0.0_dp)
        total = sum(weight)
        if (total <= 0.0_dp) error stop "normalize_diurnal_weights: curve must have positive sum"
        weight = weight / total
    end subroutine normalize_diurnal_weights

    ! Return cumulative intraday variance fractions through each bin.
    pure subroutine cumulative_diurnal_fraction(diurnal_var, cum_frac)
        real(dp), intent(in) :: diurnal_var(:)
        real(dp), intent(out) :: cum_frac(:)
        real(dp), allocatable :: weight(:)
        integer :: i

        if (size(cum_frac) /= size(diurnal_var)) then
            error stop "cumulative_diurnal_fraction: array sizes differ"
        end if
        allocate(weight(size(diurnal_var)))
        call normalize_diurnal_weights(diurnal_var, weight)
        cum_frac(1) = weight(1)
        do i = 2, size(weight)
            cum_frac(i) = cum_frac(i - 1) + weight(i)
        end do
        cum_frac(size(cum_frac)) = 1.0_dp
    end subroutine cumulative_diurnal_fraction

    ! Return the fraction of session variance from first_remaining_bin through the close.
    pure real(dp) function remaining_fraction_from_bin(diurnal_var, first_remaining_bin) result(frac)
        real(dp), intent(in) :: diurnal_var(:)
        integer, intent(in) :: first_remaining_bin
        real(dp), allocatable :: weight(:)

        allocate(weight(size(diurnal_var)))
        call normalize_diurnal_weights(diurnal_var, weight)
        if (first_remaining_bin <= 1) then
            frac = 1.0_dp
        else if (first_remaining_bin > size(weight)) then
            frac = 0.0_dp
        else
            frac = sum(weight(first_remaining_bin:))
        end if
    end function remaining_fraction_from_bin

    ! Return the fraction of session variance after last_completed_bin.
    pure real(dp) function remaining_fraction_after_bin(diurnal_var, last_completed_bin) result(frac)
        real(dp), intent(in) :: diurnal_var(:)
        integer, intent(in) :: last_completed_bin

        frac = remaining_fraction_from_bin(diurnal_var, last_completed_bin + 1)
    end function remaining_fraction_after_bin

    ! Return the clock-time session fraction from current_seconds to session_end_seconds.
    pure real(dp) function remaining_session_fraction(current_seconds, session_start_seconds, &
                                                      session_end_seconds) result(frac)
        integer, intent(in) :: current_seconds, session_start_seconds, session_end_seconds
        integer :: session_seconds, remaining_seconds

        if (session_end_seconds <= session_start_seconds) then
            error stop "remaining_session_fraction: session end must exceed start"
        end if
        session_seconds = session_end_seconds - session_start_seconds
        remaining_seconds = session_end_seconds - max(current_seconds, session_start_seconds)
        frac = real(max(remaining_seconds, 0), dp) / real(session_seconds, dp)
    end function remaining_session_fraction

    ! Return the equal-bin clock-time fraction from first_remaining_bin through the close.
    pure real(dp) function remaining_day_fraction_from_bin(nbins, first_remaining_bin) result(frac)
        integer, intent(in) :: nbins, first_remaining_bin
        integer :: remaining_bins

        if (nbins < 1) error stop "remaining_day_fraction_from_bin: nbins must be positive"
        remaining_bins = nbins - max(first_remaining_bin, 1) + 1
        frac = real(max(remaining_bins, 0), dp) / real(nbins, dp)
    end function remaining_day_fraction_from_bin

    ! Convert a full-session variance forecast to variance remaining from a diurnal curve.
    pure real(dp) function remaining_variance_from_curve(session_var, diurnal_var, &
                                                         first_remaining_bin) result(rem_var)
        real(dp), intent(in) :: session_var, diurnal_var(:)
        integer, intent(in) :: first_remaining_bin

        rem_var = max(session_var, min_var) * remaining_fraction_from_bin(diurnal_var, first_remaining_bin)
    end function remaining_variance_from_curve

    ! Return remaining variance, equal-bin day fraction, and annualized vol percent.
    pure subroutine remaining_forecast_from_curve(session_var, diurnal_var, first_remaining_bin, &
                                                  trading_days_per_year, rem_var, rem_day_frac, rem_vol_ann_pct)
        real(dp), intent(in) :: session_var, diurnal_var(:), trading_days_per_year
        integer, intent(in) :: first_remaining_bin
        real(dp), intent(out) :: rem_var, rem_day_frac, rem_vol_ann_pct

        rem_var = remaining_variance_from_curve(session_var, diurnal_var, first_remaining_bin)
        rem_day_frac = remaining_day_fraction_from_bin(size(diurnal_var), first_remaining_bin)
        rem_vol_ann_pct = trading_day_annualized_vol_pct(rem_var, rem_day_frac, trading_days_per_year)
    end subroutine remaining_forecast_from_curve

    ! Sum per-bin variance forecasts remaining in each observation's trading day.
    pure subroutine remaining_variance_from_path(h, date_id, rem_var, include_current)
        real(dp), intent(in) :: h(:)
        integer, intent(in) :: date_id(:)
        real(dp), intent(out) :: rem_var(:)
        logical, optional, intent(in) :: include_current
        integer :: i, j, first
        logical :: incl

        if (size(h) /= size(date_id) .or. size(rem_var) /= size(h)) then
            error stop "remaining_variance_from_path: array sizes differ"
        end if
        incl = .true.
        if (present(include_current)) incl = include_current
        do i = 1, size(h)
            first = i
            if (.not. incl) first = i + 1
            rem_var(i) = 0.0_dp
            do j = first, size(h)
                if (j < 1) cycle
                if (date_id(j) /= date_id(i)) exit
                rem_var(i) = rem_var(i) + max(h(j), min_var)
            end do
        end do
    end subroutine remaining_variance_from_path

    ! Sum per-bin variance forecasts over each observation's full trading day.
    pure subroutine full_session_variance_from_path(h, date_id, session_var)
        real(dp), intent(in) :: h(:)
        integer, intent(in) :: date_id(:)
        real(dp), intent(out) :: session_var(:)
        integer :: i, first, last

        if (size(h) /= size(date_id) .or. size(session_var) /= size(h)) then
            error stop "full_session_variance_from_path: array sizes differ"
        end if
        first = 1
        do while (first <= size(h))
            last = first
            do while (last < size(h) .and. date_id(last + 1) == date_id(first))
                last = last + 1
            end do
            do i = first, last
                session_var(i) = sum(max(h(first:last), min_var))
            end do
            first = last + 1
        end do
    end subroutine full_session_variance_from_path

    ! Convert variance over a horizon measured in years to annualized volatility percent.
    pure real(dp) function annualized_vol_pct(period_var, horizon_years) result(vol_pct)
        real(dp), intent(in) :: period_var, horizon_years

        if (horizon_years <= 0.0_dp) then
            vol_pct = 0.0_dp
        else
            vol_pct = 100.0_dp * sqrt(max(period_var, 0.0_dp) / horizon_years)
        end if
    end function annualized_vol_pct

    ! Convert variance over part of a trading day to trading-day annualized vol percent.
    pure real(dp) function trading_day_annualized_vol_pct(period_var, day_fraction, &
                                                          trading_days_per_year) result(vol_pct)
        real(dp), intent(in) :: period_var, day_fraction, trading_days_per_year

        if (day_fraction <= 0.0_dp .or. trading_days_per_year <= 0.0_dp) then
            vol_pct = 0.0_dp
        else
            vol_pct = annualized_vol_pct(period_var, day_fraction / trading_days_per_year)
        end if
    end function trading_day_annualized_vol_pct

    ! Combine remaining same-day variance with post-close variance to an option expiry.
    pure real(dp) function compose_expiry_variance(remaining_intraday_var, post_close_var) result(total_var)
        real(dp), intent(in) :: remaining_intraday_var, post_close_var

        total_var = max(remaining_intraday_var, 0.0_dp) + max(post_close_var, 0.0_dp)
    end function compose_expiry_variance

end module remaining_vol_mod

! Utilities for implied volatility index series.

module implied_vol_utils_mod
    use kind_mod, only: dp
    use date_mod, only: date_t, date_from_yyyymmdd, operator(-)
    implicit none
    private

    public :: adjust_implied_vol_weekday_levels

contains

    subroutine adjust_implied_vol_weekday_levels(dates, iv_raw, iv_adj)
        ! Remove multiplicative weekday seasonality from positive implied-vol levels.
        integer, intent(in) :: dates(:)
        real(dp), intent(in) :: iv_raw(:)
        real(dp), intent(out) :: iv_adj(:)
        real(dp) :: weekday_sum(5), weekday_mean(5), total_sum, total_mean
        integer :: weekday_count(5)
        integer :: i, dow, nvalid

        if (size(dates) /= size(iv_raw) .or. size(iv_raw) /= size(iv_adj)) then
            error stop "adjust_implied_vol_weekday_levels: array sizes differ"
        end if
        weekday_sum = 0.0_dp
        weekday_mean = 0.0_dp
        weekday_count = 0
        total_sum = 0.0_dp
        nvalid = 0
        do i = 1, size(iv_raw)
            if (iv_raw(i) <= 0.0_dp) cycle
            dow = day_of_week(dates(i))
            if (dow < 1 .or. dow > 5) cycle
            weekday_sum(dow) = weekday_sum(dow) + log(iv_raw(i))
            weekday_count(dow) = weekday_count(dow) + 1
            total_sum = total_sum + log(iv_raw(i))
            nvalid = nvalid + 1
        end do
        if (nvalid < 1) then
            iv_adj = iv_raw
            return
        end if
        total_mean = total_sum / real(nvalid, dp)
        do dow = 1, 5
            if (weekday_count(dow) > 0) then
                weekday_mean(dow) = weekday_sum(dow) / real(weekday_count(dow), dp)
            else
                weekday_mean(dow) = total_mean
            end if
        end do
        do i = 1, size(iv_raw)
            if (iv_raw(i) <= 0.0_dp) then
                iv_adj(i) = iv_raw(i)
            else
                dow = day_of_week(dates(i))
                if (dow >= 1 .and. dow <= 5) then
                    iv_adj(i) = exp(log(iv_raw(i)) - weekday_mean(dow) + total_mean)
                else
                    iv_adj(i) = iv_raw(i)
                end if
            end if
        end do
    end subroutine adjust_implied_vol_weekday_levels

    pure integer function day_of_week(yyyymmdd_date)
        ! Weekday number with Monday=1, ..., Sunday=7.
        integer, intent(in) :: yyyymmdd_date
        type(date_t) :: x, monday

        x = date_from_yyyymmdd(yyyymmdd_date)
        monday = date_t(1970, 1, 5)
        day_of_week = modulo(x - monday, 7) + 1
    end function day_of_week

end module implied_vol_utils_mod

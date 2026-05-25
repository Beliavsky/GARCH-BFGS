! Generic utilities for transforming and comparing dated real-valued series.

module time_series_compare_mod
    use kind_mod,    only: dp
    use strings_mod, only: uppercase
    implicit none
    private
    public :: transform_series, lag_series_one, aligned_corr, is_change_transform

contains

    subroutine transform_series(raw_dates, raw_values, transform, horizon, out_dates, out_values)
        integer, intent(in) :: raw_dates(:)
        real(dp), intent(in) :: raw_values(:)
        character(len=*), intent(in) :: transform
        integer, intent(in) :: horizon
        integer, allocatable, intent(out) :: out_dates(:)
        real(dp), allocatable, intent(out) :: out_values(:)
        integer :: n, t, h
        character(len=16) :: tr

        n = size(raw_values)
        tr = uppercase(trim(transform))
        select case (trim(tr))
        case ("LEVEL")
            allocate(out_dates(n), out_values(n))
            out_dates = raw_dates
            out_values = raw_values
        case ("LOG")
            allocate(out_dates(n), out_values(n))
            out_dates = raw_dates
            out_values = log(raw_values)
        case ("DIFF")
            h = max(horizon, 1)
            allocate(out_dates(max(n-h, 0)), out_values(max(n-h, 0)))
            do t = h + 1, n
                out_dates(t-h) = raw_dates(t)
                out_values(t-h) = raw_values(t) - raw_values(t-h)
            end do
        case ("LOG_DIFF")
            h = max(horizon, 1)
            allocate(out_dates(max(n-h, 0)), out_values(max(n-h, 0)))
            do t = h + 1, n
                out_dates(t-h) = raw_dates(t)
                out_values(t-h) = log(raw_values(t)) - log(raw_values(t-h))
            end do
        case default
            print '(A,A)', "Unknown series transform: ", trim(transform)
            error stop
        end select
    end subroutine transform_series

    subroutine lag_series_one(series_dates, series_values)
        integer, allocatable, intent(inout) :: series_dates(:)
        real(dp), allocatable, intent(inout) :: series_values(:)
        integer, allocatable :: lagged_dates(:)
        real(dp), allocatable :: lagged_values(:)
        integer :: n

        n = size(series_values)
        if (n < 2) then
            deallocate(series_dates, series_values)
            allocate(series_dates(0), series_values(0))
            return
        end if

        allocate(lagged_dates(n-1), lagged_values(n-1))
        lagged_dates = series_dates(2:n)
        lagged_values = series_values(1:n-1)
        call move_alloc(lagged_dates, series_dates)
        call move_alloc(lagged_values, series_values)
    end subroutine lag_series_one

    subroutine aligned_corr(xindex, x, yindex, y, corr, nmatch, first_index, last_index)
        integer, intent(in) :: xindex(:), yindex(:)
        real(dp), intent(in) :: x(:), y(:)
        real(dp), intent(out) :: corr
        integer, intent(out) :: nmatch, first_index, last_index
        integer :: ix, iy
        real(dp) :: sx, sy, sxx, syy, sxy, dx, dy, rn, denom

        ix = 1
        iy = 1
        nmatch = 0
        first_index = 0
        last_index = 0
        sx = 0.0_dp
        sy = 0.0_dp
        sxx = 0.0_dp
        syy = 0.0_dp
        sxy = 0.0_dp

        do while (ix <= size(xindex) .and. iy <= size(yindex))
            if (xindex(ix) == yindex(iy)) then
                nmatch = nmatch + 1
                if (nmatch == 1) first_index = xindex(ix)
                last_index = xindex(ix)
                sx = sx + x(ix)
                sy = sy + y(iy)
                sxx = sxx + x(ix)**2
                syy = syy + y(iy)**2
                sxy = sxy + x(ix)*y(iy)
                ix = ix + 1
                iy = iy + 1
            else if (xindex(ix) < yindex(iy)) then
                ix = ix + 1
            else
                iy = iy + 1
            end if
        end do

        corr = 0.0_dp
        if (nmatch <= 1) return
        rn = real(nmatch, dp)
        dx = sxx - sx*sx/rn
        dy = syy - sy*sy/rn
        denom = sqrt(max(dx, 0.0_dp) * max(dy, 0.0_dp))
        if (denom > 0.0_dp) corr = (sxy - sx*sy/rn) / denom
    end subroutine aligned_corr

    logical function is_change_transform(transform)
        character(len=*), intent(in) :: transform
        character(len=16) :: tr

        tr = uppercase(trim(transform))
        is_change_transform = trim(tr) == "DIFF" .or. trim(tr) == "LOG_DIFF"
    end function is_change_transform

end module time_series_compare_mod

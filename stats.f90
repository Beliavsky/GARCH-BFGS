! Sample statistics for real(dp) arrays.
!
! mean: sample mean
! variance: sample variance (denominator n-1)
! sd:       sample standard deviation (denominator n-1)
! sort_real: ascending in-place sort
! autocorr:  autocorrelations at lags 1:size(acf)
! print_acf_table: formatted table of ACF rows

module stats_mod
    use kind_mod, only: dp
    implicit none
    private
    public :: mean, variance, sd, sort_real, autocorr, print_acf_table

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

end module stats_mod

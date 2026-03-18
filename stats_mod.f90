! Sample statistics for real(dp) arrays.
!
! mean: sample mean
! sd:   sample standard deviation (denominator n-1)

module stats_mod
    use kind_mod, only: dp
    implicit none
    private
    public :: mean, sd

contains

    pure function mean(x) result(m)
        ! Sample mean of x.
        real(dp), intent(in) :: x(:)  ! input array
        real(dp) :: m
        m = sum(x) / size(x)
    end function mean

    pure function sd(x) result(s)
        ! Sample standard deviation of x (denominator n-1).
        real(dp), intent(in) :: x(:)  ! input array
        real(dp) :: s
        s = sqrt(sum((x - mean(x))**2) / (size(x) - 1))
    end function sd

end module stats_mod

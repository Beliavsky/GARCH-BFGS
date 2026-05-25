! Integer ranking utilities.
!
! rank_desc: rank 1 = largest element (use for log-likelihood: higher is better)
! rank_asc:  rank 1 = smallest element (use for AIC/BIC: lower is better)
!
! Ties receive the same rank (standard competition ranking).

module rank_mod
    use kind_mod, only: dp
    implicit none
    private
    public :: rank_desc, rank_asc

contains

    subroutine rank_desc(x, r)
        ! Rank x descending: r(i) = 1 means x(i) is the largest value.
        real(dp), intent(in)  :: x(:)  ! values to rank
        integer,  intent(out) :: r(:)  ! ranks (same size as x)
        integer :: i, j, n
        n = size(x)
        r = 1
        do i = 1, n
            do j = 1, n
                if (x(j) > x(i)) r(i) = r(i) + 1
            end do
        end do
    end subroutine rank_desc

    subroutine rank_asc(x, r)
        ! Rank x ascending: r(i) = 1 means x(i) is the smallest value.
        real(dp), intent(in)  :: x(:)  ! values to rank
        integer,  intent(out) :: r(:)  ! ranks (same size as x)
        integer :: i, j, n
        n = size(x)
        r = 1
        do i = 1, n
            do j = 1, n
                if (x(j) < x(i)) r(i) = r(i) + 1
            end do
        end do
    end subroutine rank_asc

end module rank_mod

! Small linear algebra module: Cholesky factorization, log-determinant, solve.
! All routines work on dense n x n matrices stored column-major (Fortran default).

module linalg_mod
    use kind_mod, only: dp
    implicit none
    private

    public :: chol_factor, chol_logdet, chol_solve_vec

contains

    ! Lower-triangular Cholesky A = L L'.  Sets ok=.false. if A is not SPD.
    subroutine chol_factor(a, n, L, ok)
        integer,  intent(in)  :: n
        real(dp), intent(in)  :: a(n,n)
        real(dp), intent(out) :: L(n,n)
        logical,  intent(out) :: ok
        integer  :: i, j
        real(dp) :: s
        L  = 0.0_dp
        ok = .true.
        do j = 1, n
            s = a(j,j) - dot_product(L(j,1:j-1), L(j,1:j-1))
            if (s <= 0.0_dp) then
                ok = .false.
                return
            end if
            L(j,j) = sqrt(s)
            do i = j+1, n
                L(i,j) = (a(i,j) - dot_product(L(i,1:j-1), L(j,1:j-1))) / L(j,j)
            end do
        end do
    end subroutine chol_factor

    ! log|A| = 2 * sum_i log(L_ii) from the Cholesky factor L.
    function chol_logdet(L, n) result(ld)
        integer,  intent(in) :: n
        real(dp), intent(in) :: L(n,n)
        real(dp) :: ld
        integer  :: i
        ld = 0.0_dp
        do i = 1, n
            ld = ld + log(L(i,i))
        end do
        ld = 2.0_dp * ld
    end function chol_logdet

    ! Solve A x = b where A = L L' (L lower-triangular Cholesky factor).
    subroutine chol_solve_vec(L, n, b, x)
        integer,  intent(in)  :: n
        real(dp), intent(in)  :: L(n,n), b(n)
        real(dp), intent(out) :: x(n)
        real(dp) :: y(n)
        integer  :: i
        ! Forward: L y = b
        do i = 1, n
            y(i) = (b(i) - dot_product(L(i,1:i-1), y(1:i-1))) / L(i,i)
        end do
        ! Backward: L' x = y
        do i = n, 1, -1
            x(i) = (y(i) - dot_product(L(i+1:n,i), x(i+1:n))) / L(i,i)
        end do
    end subroutine chol_solve_vec

end module linalg_mod

! BFGS quasi-Newton minimiser with backtracking Armijo line search.
!
! The inverse Hessian approximation H is initialised to the identity and
! updated via the standard BFGS formula.  When the curvature condition
! s'y > 0 fails the update is skipped.  When the search direction is not
! a descent direction H is reset to the identity.

module bfgs_mod
    use kind_mod, only: dp
    implicit none
    private

    public :: bfgs_minimize

    abstract interface
        subroutine obj_iface(p, np, f, g)
            import :: dp
            integer,  intent(in)  :: np
            real(dp), intent(in)  :: p(np)
            real(dp), intent(out) :: f, g(np)
        end subroutine obj_iface
    end interface

contains

    subroutine bfgs_minimize(func, p, np, max_iter, gtol, f_opt, n_iter, converged)
        procedure(obj_iface)    :: func
        integer,  intent(in)    :: np, max_iter
        real(dp), intent(inout) :: p(np)
        real(dp), intent(in)    :: gtol
        real(dp), intent(out)   :: f_opt
        integer,  intent(out)   :: n_iter
        logical,  intent(out)   :: converged

        real(dp) :: f, f_new
        real(dp), allocatable :: g(:), p_new(:), g_new(:)
        real(dp), allocatable :: H(:,:), d(:), s(:), y(:)
        real(dp) :: sy, rho
        real(dp), allocatable :: A(:,:), tmp(:,:)
        integer  :: i, j, iter

        allocate(g(np), p_new(np), g_new(np), H(np,np), d(np), s(np), y(np), A(np,np), tmp(np,np))

        ! initialise inverse Hessian approximation to identity
        H = 0.0_dp
        do i = 1, np
            H(i,i) = 1.0_dp
        end do

        call func(p, np, f, g)

        converged = .false.
        n_iter    = 0

        do iter = 1, max_iter
            n_iter = iter

            if (norm2(g) < gtol) then
                converged = .true.
                exit
            end if

            ! search direction d = -H * g
            do i = 1, np
                d(i) = -sum(H(i,:) * g)
            end do

            ! if d is not a descent direction, reset H and use steepest descent
            if (dot_product(g, d) >= 0.0_dp) then
                H = 0.0_dp
                do i = 1, np
                    H(i,i) = 1.0_dp
                end do
                d = -g
            end if

            ! backtracking line search (Armijo condition)
            call backtrack(func, p, f, g, d, np, p_new, f_new, g_new)

            s  = p_new - p
            y  = g_new - g
            sy = dot_product(s, y)

            p = p_new
            f = f_new
            g = g_new

            ! BFGS inverse Hessian update
            ! H <- (I - rho*s*y') * H * (I - rho*y*s') + rho*s*s'
            if (sy > 1.0e-10_dp * norm2(s) * norm2(y)) then
                rho = 1.0_dp / sy

                ! A = I - rho * outer(s, y)
                do i = 1, np
                    do j = 1, np
                        A(i,j) = -rho * s(i) * y(j)
                    end do
                    A(i,i) = A(i,i) + 1.0_dp
                end do

                ! tmp = H * A^T
                do i = 1, np
                    do j = 1, np
                        tmp(i,j) = sum(H(i,:) * A(j,:))
                    end do
                end do

                ! H = A * tmp + rho * outer(s, s)
                do i = 1, np
                    do j = 1, np
                        H(i,j) = sum(A(i,:) * tmp(:,j)) + rho * s(i) * s(j)
                    end do
                end do
            end if

        end do

        f_opt = f
        deallocate(g, p_new, g_new, H, d, s, y, A, tmp)
    end subroutine bfgs_minimize

    ! Backtracking line search: halve alpha until the Armijo condition holds
    ! or alpha falls below a floor (in which case we accept and move on).
    subroutine backtrack(func, p, f, g, d, np, p_new, f_new, g_new)
        procedure(obj_iface)  :: func
        integer,  intent(in)  :: np
        real(dp), intent(in)  :: p(np), f, g(np), d(np)
        real(dp), intent(out) :: p_new(np), f_new, g_new(np)

        real(dp), parameter :: c1        = 1.0e-4_dp
        real(dp), parameter :: rho_ls    = 0.5_dp
        real(dp), parameter :: alpha_min = 1.0e-12_dp
        integer,  parameter :: max_ls    = 60

        real(dp) :: alpha, slope
        integer  :: ls_iter

        slope = dot_product(g, d)
        alpha = 1.0_dp

        do ls_iter = 1, max_ls
            p_new = p + alpha * d
            call func(p_new, np, f_new, g_new)
            if (f_new <= f + c1 * alpha * slope) return
            if (alpha < alpha_min) return
            alpha = rho_ls * alpha
        end do
    end subroutine backtrack

end module bfgs_mod

! Multiplicative component sGARCH for intraday returns.
!
! The model follows the rugarch mcsGARCH structure:
!   r_t = sqrt(daily_var_t * diurnal_var_t * q_t) * z_t
!   q_t = omega + alpha * e_{t-1}^2 + beta * q_{t-1}
! where e_t = r_t / sqrt(daily_var_t * diurnal_var_t).

module garch_mcsgarch_mod
    use kind_mod, only: dp
    use bfgs_mod, only: bfgs_minimize
    use stats_mod, only: variance
    implicit none
    private

    real(dp), parameter :: min_var = 1.0e-12_dp
    real(dp), parameter :: min_pdf = 1.0e-300_dp
    real(dp), parameter :: persist_max = 0.999_dp
    real(dp), parameter :: log_sqrt_2pi = 0.91893853320467274178_dp

    type, public :: mcsgarch_params_t
        real(dp) :: omega = 0.02_dp
        real(dp) :: alpha = 0.06_dp
        real(dp) :: beta = 0.90_dp
    end type mcsgarch_params_t

    type, public :: mcsgarch_fit_result_t
        type(mcsgarch_params_t) :: params
        real(dp) :: nll = huge(1.0_dp)
        real(dp) :: loglik = -huge(1.0_dp)
        real(dp) :: persist = 0.0_dp
        integer :: niter = 0
        logical :: converged = .false.
    end type mcsgarch_fit_result_t

    real(dp), allocatable, save :: obj_y(:), obj_daily_var(:), obj_diurnal_var(:)

    public :: fit_mcsgarch
    public :: mcsgarch_filter, mcsgarch_simulate
    public :: estimate_diurnal_variance

contains

    subroutine fit_mcsgarch(y, daily_var, bin_id, max_iter, gtol, result, diurnal_var, q)
        real(dp), intent(in) :: y(:), daily_var(:), gtol
        integer, intent(in) :: bin_id(:), max_iter
        type(mcsgarch_fit_result_t), intent(out) :: result
        real(dp), optional, intent(out) :: diurnal_var(:), q(:)
        real(dp), allocatable :: p(:), grad(:), dvar(:), qtmp(:)
        real(dp) :: fbest, vnorm
        integer :: niter
        logical :: converged

        call check_inputs(y, daily_var, bin_id)
        allocate(dvar(size(y)), qtmp(size(y)), p(3), grad(3))
        call estimate_diurnal_variance(y, daily_var, bin_id, dvar)

        if (allocated(obj_y)) deallocate(obj_y)
        if (allocated(obj_daily_var)) deallocate(obj_daily_var)
        if (allocated(obj_diurnal_var)) deallocate(obj_diurnal_var)
        allocate(obj_y(size(y)), obj_daily_var(size(y)), obj_diurnal_var(size(y)))
        obj_y = y
        obj_daily_var = max(daily_var, min_var)
        obj_diurnal_var = max(dvar, min_var)

        vnorm = max(variance(y / sqrt(obj_daily_var * obj_diurnal_var)), min_var)
        call pack_params(mcsgarch_params_t(vnorm * 0.04_dp, 0.06_dp, 0.90_dp), p)
        call bfgs_minimize(mcsgarch_obj, p, 3, max_iter, gtol, fbest, niter, converged)
        call unpack_params(p, result%params)
        call mcsgarch_obj(p, 3, result%nll, grad)
        call mcsgarch_filter(y, daily_var, dvar, result%params, qtmp)

        result%loglik = -result%nll * real(size(y), dp)
        result%persist = result%params%alpha + result%params%beta
        result%niter = niter
        result%converged = converged
        if (present(diurnal_var)) diurnal_var = dvar
        if (present(q)) q = qtmp

        deallocate(p, grad, dvar, qtmp)
    end subroutine fit_mcsgarch

    subroutine estimate_diurnal_variance(y, daily_var, bin_id, diurnal_var)
        real(dp), intent(in) :: y(:), daily_var(:)
        integer, intent(in) :: bin_id(:)
        real(dp), intent(out) :: diurnal_var(:)
        real(dp), allocatable :: bin_sum(:)
        integer, allocatable :: bin_count(:)
        integer :: nbins, i, b
        real(dp) :: avg

        call check_inputs(y, daily_var, bin_id)
        nbins = maxval(bin_id)
        if (minval(bin_id) < 1) error stop "estimate_diurnal_variance: bin_id must be positive"
        allocate(bin_sum(nbins), bin_count(nbins))
        bin_sum = 0.0_dp
        bin_count = 0
        do i = 1, size(y)
            b = bin_id(i)
            bin_sum(b) = bin_sum(b) + y(i)**2 / max(daily_var(i), min_var)
            bin_count(b) = bin_count(b) + 1
        end do
        avg = 0.0_dp
        do b = 1, nbins
            if (bin_count(b) > 0) then
                bin_sum(b) = bin_sum(b) / real(bin_count(b), dp)
            else
                bin_sum(b) = 1.0_dp
            end if
            avg = avg + bin_sum(b)
        end do
        avg = max(avg / real(nbins, dp), min_var)
        bin_sum = max(bin_sum / avg, min_var)
        do i = 1, size(y)
            diurnal_var(i) = bin_sum(bin_id(i))
        end do
        deallocate(bin_sum, bin_count)
    end subroutine estimate_diurnal_variance

    subroutine mcsgarch_filter(y, daily_var, diurnal_var, params, q)
        real(dp), intent(in) :: y(:), daily_var(:), diurnal_var(:)
        type(mcsgarch_params_t), intent(in) :: params
        real(dp), intent(out) :: q(:)
        real(dp), allocatable :: e(:)
        integer :: t

        if (size(q) /= size(y)) error stop "mcsgarch_filter: q has wrong size"
        if (.not. params_valid(params)) error stop "mcsgarch_filter: invalid parameters"
        allocate(e(size(y)))
        e = y / sqrt(max(daily_var, min_var) * max(diurnal_var, min_var))
        q(1) = max(params%omega / max(1.0_dp - params%alpha - params%beta, 1.0e-8_dp), variance(e))
        q(1) = max(q(1), min_var)
        do t = 2, size(y)
            q(t) = max(params%omega + params%alpha*e(t - 1)**2 + params%beta*q(t - 1), min_var)
        end do
        deallocate(e)
    end subroutine mcsgarch_filter

    subroutine mcsgarch_simulate(params, daily_var, diurnal_var, seed_val, y, q)
        type(mcsgarch_params_t), intent(in) :: params
        real(dp), intent(in) :: daily_var(:), diurnal_var(:)
        integer, intent(in) :: seed_val
        real(dp), intent(out) :: y(:)
        real(dp), optional, intent(out) :: q(:)
        real(dp) :: qprev, eprev, qnow, z
        integer :: t

        if (size(y) /= size(daily_var) .or. size(y) /= size(diurnal_var)) then
            error stop "mcsgarch_simulate: array sizes differ"
        end if
        if (.not. params_valid(params)) error stop "mcsgarch_simulate: invalid parameters"
        call seed_rng(seed_val)
        qprev = max(params%omega / max(1.0_dp - params%alpha - params%beta, 1.0e-8_dp), min_var)
        eprev = 0.0_dp
        do t = 1, size(y)
            if (t == 1) then
                qnow = qprev
            else
                qnow = max(params%omega + params%alpha*eprev**2 + params%beta*qprev, min_var)
            end if
            z = random_normal()
            eprev = sqrt(qnow) * z
            y(t) = eprev * sqrt(max(daily_var(t), min_var) * max(diurnal_var(t), min_var))
            if (present(q)) q(t) = qnow
            qprev = qnow
        end do
    end subroutine mcsgarch_simulate

    subroutine mcsgarch_obj(p, np, f, g)
        integer, intent(in) :: np
        real(dp), intent(in) :: p(np)
        real(dp), intent(out) :: f, g(np)
        real(dp), allocatable :: pp(:), pm(:)
        real(dp) :: step, fp, fm
        integer :: j

        f = mcsgarch_nll_from_p(p)
        allocate(pp(np), pm(np))
        do j = 1, np
            step = 1.0e-5_dp * max(1.0_dp, abs(p(j)))
            pp = p
            pm = p
            pp(j) = pp(j) + step
            pm(j) = pm(j) - step
            fp = mcsgarch_nll_from_p(pp)
            fm = mcsgarch_nll_from_p(pm)
            g(j) = (fp - fm) / (2.0_dp * step)
        end do
        deallocate(pp, pm)
    end subroutine mcsgarch_obj

    real(dp) function mcsgarch_nll_from_p(p)
        real(dp), intent(in) :: p(:)
        type(mcsgarch_params_t) :: params
        real(dp), allocatable :: q(:)
        real(dp) :: h, z, loss
        integer :: t

        call unpack_params(p, params)
        if (.not. params_valid(params)) then
            mcsgarch_nll_from_p = huge(1.0_dp) / 10.0_dp
            return
        end if
        allocate(q(size(obj_y)))
        call mcsgarch_filter(obj_y, obj_daily_var, obj_diurnal_var, params, q)
        loss = 0.0_dp
        do t = 1, size(obj_y)
            h = max(obj_daily_var(t) * obj_diurnal_var(t) * q(t), min_var)
            z = obj_y(t) / sqrt(h)
            loss = loss + log_sqrt_2pi + 0.5_dp*log(h) + 0.5_dp*z*z
        end do
        mcsgarch_nll_from_p = loss / real(size(obj_y), dp)
        deallocate(q)
    end function mcsgarch_nll_from_p

    subroutine pack_params(params, p)
        type(mcsgarch_params_t), intent(in) :: params
        real(dp), intent(out) :: p(3)
        real(dp) :: rest, a, b

        p(1) = log(max(params%omega, min_var))
        a = max(params%alpha, min_var)
        b = max(params%beta, min_var)
        rest = max(persist_max - a - b, min_var)
        p(2) = log(a / rest)
        p(3) = log(b / rest)
    end subroutine pack_params

    subroutine unpack_params(p, params)
        real(dp), intent(in) :: p(:)
        type(mcsgarch_params_t), intent(out) :: params
        real(dp) :: ea, eb, den

        params%omega = exp(max(min(p(1), 50.0_dp), -50.0_dp))
        ea = exp(max(min(p(2), 50.0_dp), -50.0_dp))
        eb = exp(max(min(p(3), 50.0_dp), -50.0_dp))
        den = 1.0_dp + ea + eb
        params%alpha = persist_max * ea / den
        params%beta = persist_max * eb / den
    end subroutine unpack_params

    logical function params_valid(params)
        type(mcsgarch_params_t), intent(in) :: params

        params_valid = params%omega > 0.0_dp .and. params%alpha >= 0.0_dp .and. &
                       params%beta >= 0.0_dp .and. params%alpha + params%beta < 1.0_dp
    end function params_valid

    subroutine check_inputs(y, daily_var, bin_id)
        real(dp), intent(in) :: y(:), daily_var(:)
        integer, intent(in) :: bin_id(:)

        if (size(y) < 3) error stop "mcsGARCH: need at least three observations"
        if (size(daily_var) /= size(y) .or. size(bin_id) /= size(y)) then
            error stop "mcsGARCH: y, daily_var, and bin_id sizes differ"
        end if
        if (any(daily_var <= 0.0_dp)) error stop "mcsGARCH: daily_var must be positive"
    end subroutine check_inputs

    subroutine seed_rng(seed_val)
        integer, intent(in) :: seed_val
        integer :: sz
        integer, allocatable :: seed_arr(:)

        call random_seed(size=sz)
        allocate(seed_arr(sz))
        seed_arr = seed_val
        call random_seed(put=seed_arr)
        deallocate(seed_arr)
    end subroutine seed_rng

    real(dp) function random_normal()
        real(dp) :: u1, u2
        real(dp), parameter :: two_pi = 2.0_dp*acos(-1.0_dp)

        do
            call random_number(u1)
            if (u1 > 0.0_dp) exit
        end do
        call random_number(u2)
        random_normal = sqrt(-2.0_dp * log(u1)) * cos(two_pi*u2)
    end function random_normal

end module garch_mcsgarch_mod

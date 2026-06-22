! Multiplicative component sGARCH for intraday returns.
!
! The model follows the rugarch mcsGARCH structure:
!   r_t = sqrt(daily_var_t * diurnal_var_t * q_t) * z_t
!   q_t = omega + alpha * e_{t-1}^2 + beta * q_{t-1}
! where e_t = r_t / sqrt(daily_var_t * diurnal_var_t).
! NAGARCH-style and GJR-style asymmetric q_t recursions are also supported.

module garch_mcsgarch_mod
    use kind_mod, only: dp
    use math_const_mod, only: log_sqrt_2pi, pi
    use distributions_mod, only: pdf_fs_skewt
    use bfgs_mod, only: bfgs_minimize
    use stats_mod, only: variance
    implicit none
    private

    real(dp), parameter :: min_var = 1.0e-12_dp
    real(dp), parameter :: min_pdf = 1.0e-300_dp
    real(dp), parameter :: persist_max = 0.999_dp
    integer, parameter :: model_sym = 1
    integer, parameter :: model_nagarch = 2
    integer, parameter :: model_gjr = 3
    integer, parameter :: dist_normal = 1
    integer, parameter :: dist_t = 2
    integer, parameter :: dist_fs_skewt = 3

    type, public :: mcsgarch_params_t
        real(dp) :: omega = 0.02_dp
        real(dp) :: alpha = 0.06_dp
        real(dp) :: beta = 0.90_dp
        real(dp) :: theta = 0.0_dp
        real(dp) :: gamma = 0.0_dp
    end type mcsgarch_params_t

    type, public :: mcsgarch_fit_result_t
        type(mcsgarch_params_t) :: params
        real(dp) :: nll = huge(1.0_dp)
        real(dp) :: loglik = -huge(1.0_dp)
        real(dp) :: persist = 0.0_dp
        real(dp) :: nu = 0.0_dp
        real(dp) :: xi = 0.0_dp
        integer :: niter = 0
        logical :: converged = .false.
    end type mcsgarch_fit_result_t

    real(dp), allocatable, save :: obj_y(:), obj_daily_var(:), obj_diurnal_var(:)
    integer, save :: obj_model = model_sym
    integer, save :: obj_dist = dist_normal

    public :: fit_mcsgarch
    public :: fit_mcsgarch_nagarch
    public :: fit_mcsgarch_gjr
    public :: fit_mcsgarch_t
    public :: fit_mcsgarch_nagarch_t
    public :: fit_mcsgarch_gjr_t
    public :: fit_mcsgarch_fs_skewt
    public :: fit_mcsgarch_nagarch_fs_skewt
    public :: fit_mcsgarch_gjr_fs_skewt
    public :: mcsgarch_filter, mcsgarch_simulate
    public :: mcsgarch_nagarch_filter
    public :: mcsgarch_gjr_filter
    public :: estimate_diurnal_variance
    public :: mcsgarch_remaining_var

contains

    ! Fit symmetric MCS-GARCH intraday dynamics.
    subroutine fit_mcsgarch(y, daily_var, bin_id, max_iter, gtol, result, diurnal_var, q, &
                            smooth_diurnal, smooth_half_width)
        real(dp), intent(in) :: y(:), daily_var(:), gtol
        integer, intent(in) :: bin_id(:), max_iter
        type(mcsgarch_fit_result_t), intent(out) :: result
        real(dp), optional, intent(out) :: diurnal_var(:), q(:)
        logical, optional, intent(in) :: smooth_diurnal
        integer, optional, intent(in) :: smooth_half_width

        call fit_mcsgarch_core(y, daily_var, bin_id, max_iter, gtol, model_sym, dist_normal, result, diurnal_var, q, &
                               smooth_diurnal, smooth_half_width)
    end subroutine fit_mcsgarch

    ! Fit NAGARCH-style asymmetric MCS-GARCH intraday dynamics.
    subroutine fit_mcsgarch_nagarch(y, daily_var, bin_id, max_iter, gtol, result, diurnal_var, q, &
                                    smooth_diurnal, smooth_half_width)
        real(dp), intent(in) :: y(:), daily_var(:), gtol
        integer, intent(in) :: bin_id(:), max_iter
        type(mcsgarch_fit_result_t), intent(out) :: result
        real(dp), optional, intent(out) :: diurnal_var(:), q(:)
        logical, optional, intent(in) :: smooth_diurnal
        integer, optional, intent(in) :: smooth_half_width

        call fit_mcsgarch_core(y, daily_var, bin_id, max_iter, gtol, model_nagarch, dist_normal, result, diurnal_var, q, &
                               smooth_diurnal, smooth_half_width)
    end subroutine fit_mcsgarch_nagarch

    ! Fit GJR-style asymmetric MCS-GARCH intraday dynamics.
    subroutine fit_mcsgarch_gjr(y, daily_var, bin_id, max_iter, gtol, result, diurnal_var, q, &
                                smooth_diurnal, smooth_half_width)
        real(dp), intent(in) :: y(:), daily_var(:), gtol
        integer, intent(in) :: bin_id(:), max_iter
        type(mcsgarch_fit_result_t), intent(out) :: result
        real(dp), optional, intent(out) :: diurnal_var(:), q(:)
        logical, optional, intent(in) :: smooth_diurnal
        integer, optional, intent(in) :: smooth_half_width

        call fit_mcsgarch_core(y, daily_var, bin_id, max_iter, gtol, model_gjr, dist_normal, result, diurnal_var, q, &
                               smooth_diurnal, smooth_half_width)
    end subroutine fit_mcsgarch_gjr

    ! Fit symmetric MCS-GARCH with Student-t innovations.
    subroutine fit_mcsgarch_t(y, daily_var, bin_id, max_iter, gtol, result, diurnal_var, q, &
                              smooth_diurnal, smooth_half_width)
        real(dp), intent(in) :: y(:), daily_var(:), gtol
        integer, intent(in) :: bin_id(:), max_iter
        type(mcsgarch_fit_result_t), intent(out) :: result
        real(dp), optional, intent(out) :: diurnal_var(:), q(:)
        logical, optional, intent(in) :: smooth_diurnal
        integer, optional, intent(in) :: smooth_half_width

        call fit_mcsgarch_core(y, daily_var, bin_id, max_iter, gtol, model_sym, dist_t, result, diurnal_var, q, &
                               smooth_diurnal, smooth_half_width)
    end subroutine fit_mcsgarch_t

    ! Fit NAGARCH-style MCS-GARCH with Student-t innovations.
    subroutine fit_mcsgarch_nagarch_t(y, daily_var, bin_id, max_iter, gtol, result, diurnal_var, q, &
                                      smooth_diurnal, smooth_half_width)
        real(dp), intent(in) :: y(:), daily_var(:), gtol
        integer, intent(in) :: bin_id(:), max_iter
        type(mcsgarch_fit_result_t), intent(out) :: result
        real(dp), optional, intent(out) :: diurnal_var(:), q(:)
        logical, optional, intent(in) :: smooth_diurnal
        integer, optional, intent(in) :: smooth_half_width

        call fit_mcsgarch_core(y, daily_var, bin_id, max_iter, gtol, model_nagarch, dist_t, result, diurnal_var, q, &
                               smooth_diurnal, smooth_half_width)
    end subroutine fit_mcsgarch_nagarch_t

    ! Fit GJR-style MCS-GARCH with Student-t innovations.
    subroutine fit_mcsgarch_gjr_t(y, daily_var, bin_id, max_iter, gtol, result, diurnal_var, q, &
                                  smooth_diurnal, smooth_half_width)
        real(dp), intent(in) :: y(:), daily_var(:), gtol
        integer, intent(in) :: bin_id(:), max_iter
        type(mcsgarch_fit_result_t), intent(out) :: result
        real(dp), optional, intent(out) :: diurnal_var(:), q(:)
        logical, optional, intent(in) :: smooth_diurnal
        integer, optional, intent(in) :: smooth_half_width

        call fit_mcsgarch_core(y, daily_var, bin_id, max_iter, gtol, model_gjr, dist_t, result, diurnal_var, q, &
                               smooth_diurnal, smooth_half_width)
    end subroutine fit_mcsgarch_gjr_t

    ! Fit symmetric MCS-GARCH with Fernandez-Steel skewed-t innovations.
    subroutine fit_mcsgarch_fs_skewt(y, daily_var, bin_id, max_iter, gtol, result, diurnal_var, q, &
                                     smooth_diurnal, smooth_half_width)
        real(dp), intent(in) :: y(:), daily_var(:), gtol
        integer, intent(in) :: bin_id(:), max_iter
        type(mcsgarch_fit_result_t), intent(out) :: result
        real(dp), optional, intent(out) :: diurnal_var(:), q(:)
        logical, optional, intent(in) :: smooth_diurnal
        integer, optional, intent(in) :: smooth_half_width

        call fit_mcsgarch_core(y, daily_var, bin_id, max_iter, gtol, model_sym, dist_fs_skewt, result, diurnal_var, q, &
                               smooth_diurnal, smooth_half_width)
    end subroutine fit_mcsgarch_fs_skewt

    ! Fit NAGARCH-style MCS-GARCH with Fernandez-Steel skewed-t innovations.
    subroutine fit_mcsgarch_nagarch_fs_skewt(y, daily_var, bin_id, max_iter, gtol, result, diurnal_var, q, &
                                             smooth_diurnal, smooth_half_width)
        real(dp), intent(in) :: y(:), daily_var(:), gtol
        integer, intent(in) :: bin_id(:), max_iter
        type(mcsgarch_fit_result_t), intent(out) :: result
        real(dp), optional, intent(out) :: diurnal_var(:), q(:)
        logical, optional, intent(in) :: smooth_diurnal
        integer, optional, intent(in) :: smooth_half_width

        call fit_mcsgarch_core(y, daily_var, bin_id, max_iter, gtol, model_nagarch, dist_fs_skewt, result, diurnal_var, q, &
                               smooth_diurnal, smooth_half_width)
    end subroutine fit_mcsgarch_nagarch_fs_skewt

    ! Fit GJR-style MCS-GARCH with Fernandez-Steel skewed-t innovations.
    subroutine fit_mcsgarch_gjr_fs_skewt(y, daily_var, bin_id, max_iter, gtol, result, diurnal_var, q, &
                                         smooth_diurnal, smooth_half_width)
        real(dp), intent(in) :: y(:), daily_var(:), gtol
        integer, intent(in) :: bin_id(:), max_iter
        type(mcsgarch_fit_result_t), intent(out) :: result
        real(dp), optional, intent(out) :: diurnal_var(:), q(:)
        logical, optional, intent(in) :: smooth_diurnal
        integer, optional, intent(in) :: smooth_half_width

        call fit_mcsgarch_core(y, daily_var, bin_id, max_iter, gtol, model_gjr, dist_fs_skewt, result, diurnal_var, q, &
                               smooth_diurnal, smooth_half_width)
    end subroutine fit_mcsgarch_gjr_fs_skewt

    ! Shared fitter for symmetric and asymmetric intraday q_t recursions.
    subroutine fit_mcsgarch_core(y, daily_var, bin_id, max_iter, gtol, model, dist, result, diurnal_var, q, &
                                 smooth_diurnal, smooth_half_width)
        real(dp), intent(in) :: y(:), daily_var(:), gtol
        integer, intent(in) :: bin_id(:), max_iter
        integer, intent(in) :: model, dist
        type(mcsgarch_fit_result_t), intent(out) :: result
        real(dp), optional, intent(out) :: diurnal_var(:), q(:)
        logical, optional, intent(in) :: smooth_diurnal
        integer, optional, intent(in) :: smooth_half_width
        real(dp), allocatable :: p(:), grad(:), dvar(:), qtmp(:)
        real(dp) :: fbest, vnorm
        integer :: niter, np_model, np
        logical :: converged

        call check_inputs(y, daily_var, bin_id)
        if (model == model_sym) then
            np_model = 3
        else
            np_model = 4
        end if
        np = np_model
        if (dist == dist_t) np = np + 1
        if (dist == dist_fs_skewt) np = np + 2
        allocate(dvar(size(y)), qtmp(size(y)), p(np), grad(np))
        call estimate_diurnal_variance(y, daily_var, bin_id, dvar, smooth_diurnal, smooth_half_width)

        if (allocated(obj_y)) deallocate(obj_y)
        if (allocated(obj_daily_var)) deallocate(obj_daily_var)
        if (allocated(obj_diurnal_var)) deallocate(obj_diurnal_var)
        allocate(obj_y(size(y)), obj_daily_var(size(y)), obj_diurnal_var(size(y)))
        obj_y = y
        obj_daily_var = max(daily_var, min_var)
        obj_diurnal_var = max(dvar, min_var)
        obj_model = model
        obj_dist = dist

        vnorm = max(variance(y / sqrt(obj_daily_var * obj_diurnal_var)), min_var)
        call pack_params(mcsgarch_params_t(omega=vnorm * 0.04_dp, alpha=0.06_dp, beta=0.90_dp, theta=0.0_dp, &
                                           gamma=0.05_dp), model, dist, p)
        call bfgs_minimize(mcsgarch_obj, p, np, max_iter, gtol, fbest, niter, converged)
        call unpack_params(p, model, dist, result%params, result%nu, result%xi)
        call mcsgarch_obj(p, np, result%nll, grad)
        select case (model)
        case (model_nagarch)
            call mcsgarch_nagarch_filter(y, daily_var, dvar, result%params, qtmp)
        case (model_gjr)
            call mcsgarch_gjr_filter(y, daily_var, dvar, result%params, qtmp)
        case default
            call mcsgarch_filter(y, daily_var, dvar, result%params, qtmp)
        end select

        result%loglik = -result%nll * real(size(y), dp)
        result%persist = mcsgarch_persist(result%params, model)
        if (dist == dist_normal) result%nu = 0.0_dp
        if (dist /= dist_fs_skewt) result%xi = 0.0_dp
        result%niter = niter
        result%converged = converged
        if (present(diurnal_var)) diurnal_var = dvar
        if (present(q)) q = qtmp

        deallocate(p, grad, dvar, qtmp)
    end subroutine fit_mcsgarch_core

    subroutine estimate_diurnal_variance(y, daily_var, bin_id, diurnal_var, smooth_diurnal, smooth_half_width)
        real(dp), intent(in) :: y(:), daily_var(:)
        integer, intent(in) :: bin_id(:)
        real(dp), intent(out) :: diurnal_var(:)
        logical, optional, intent(in) :: smooth_diurnal
        integer, optional, intent(in) :: smooth_half_width
        real(dp), allocatable :: bin_sum(:), bin_curve(:)
        integer, allocatable :: bin_count(:)
        integer :: nbins, i, b
        real(dp) :: avg
        logical :: do_smooth
        integer :: half_width

        call check_inputs(y, daily_var, bin_id)
        nbins = maxval(bin_id)
        if (minval(bin_id) < 1) error stop "estimate_diurnal_variance: bin_id must be positive"
        allocate(bin_sum(nbins), bin_curve(nbins), bin_count(nbins))
        bin_sum = 0.0_dp
        bin_count = 0
        do_smooth = .false.
        half_width = 2
        if (present(smooth_diurnal)) do_smooth = smooth_diurnal
        if (present(smooth_half_width)) half_width = smooth_half_width
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
        bin_curve = max(bin_sum / avg, min_var)
        if (do_smooth) call smooth_diurnal_curve(bin_curve, bin_count, half_width)
        do i = 1, size(y)
            diurnal_var(i) = bin_curve(bin_id(i))
        end do
        deallocate(bin_sum, bin_curve, bin_count)
    end subroutine estimate_diurnal_variance

    subroutine smooth_diurnal_curve(bin_curve, bin_count, half_width)
        real(dp), intent(inout) :: bin_curve(:)
        integer, intent(in) :: bin_count(:), half_width
        real(dp), allocatable :: smooth(:)
        real(dp) :: wsum, lsum, weight, avg
        integer :: b, j, lo, hi, dist

        if (half_width < 1) return
        allocate(smooth(size(bin_curve)))
        smooth = bin_curve
        do b = 1, size(bin_curve)
            if (bin_count(b) < 1) cycle
            lo = max(1, b - half_width)
            hi = min(size(bin_curve), b + half_width)
            wsum = 0.0_dp
            lsum = 0.0_dp
            do j = lo, hi
                if (bin_count(j) < 1) cycle
                dist = abs(j - b)
                weight = real(half_width + 1 - dist, dp)
                wsum = wsum + weight
                lsum = lsum + weight * log(max(bin_curve(j), min_var))
            end do
            if (wsum > 0.0_dp) smooth(b) = exp(lsum / wsum)
        end do

        avg = 0.0_dp
        wsum = 0.0_dp
        do b = 1, size(smooth)
            if (bin_count(b) < 1) cycle
            avg = avg + smooth(b)
            wsum = wsum + 1.0_dp
        end do
        if (wsum > 0.0_dp) then
            avg = max(avg / wsum, min_var)
            smooth = max(smooth / avg, min_var)
        end if
        bin_curve = smooth
        deallocate(smooth)
    end subroutine smooth_diurnal_curve

    subroutine mcsgarch_filter(y, daily_var, diurnal_var, params, q)
        real(dp), intent(in) :: y(:), daily_var(:), diurnal_var(:)
        type(mcsgarch_params_t), intent(in) :: params
        real(dp), intent(out) :: q(:)
        real(dp), allocatable :: e(:)
        integer :: t

        if (size(q) /= size(y)) error stop "mcsgarch_filter: q has wrong size"
        if (.not. params_valid(params, model_sym)) error stop "mcsgarch_filter: invalid parameters"
        allocate(e(size(y)))
        e = y / sqrt(max(daily_var, min_var) * max(diurnal_var, min_var))
        q(1) = max(params%omega / max(1.0_dp - params%alpha - params%beta, 1.0e-8_dp), variance(e))
        q(1) = max(q(1), min_var)
        do t = 2, size(y)
            q(t) = max(params%omega + params%alpha*e(t - 1)**2 + params%beta*q(t - 1), min_var)
        end do
        deallocate(e)
    end subroutine mcsgarch_filter

    ! Filter q_t using a NAGARCH-style shifted news impact curve.
    subroutine mcsgarch_nagarch_filter(y, daily_var, diurnal_var, params, q)
        real(dp), intent(in) :: y(:), daily_var(:), diurnal_var(:)
        type(mcsgarch_params_t), intent(in) :: params
        real(dp), intent(out) :: q(:)
        real(dp), allocatable :: e(:)
        integer :: t

        if (size(q) /= size(y)) error stop "mcsgarch_nagarch_filter: q has wrong size"
        if (.not. params_valid(params, model_nagarch)) error stop "mcsgarch_nagarch_filter: invalid parameters"
        allocate(e(size(y)))
        e = y / sqrt(max(daily_var, min_var) * max(diurnal_var, min_var))
        q(1) = max(params%omega / max(1.0_dp - mcsgarch_persist(params, model_nagarch), 1.0e-8_dp), variance(e))
        q(1) = max(q(1), min_var)
        do t = 2, size(y)
            q(t) = max(params%omega + params%alpha*(e(t - 1) - params%theta*sqrt(q(t - 1)))**2 + &
                       params%beta*q(t - 1), min_var)
        end do
        deallocate(e)
    end subroutine mcsgarch_nagarch_filter

    ! Filter q_t using a GJR-style negative-shock news impact term.
    subroutine mcsgarch_gjr_filter(y, daily_var, diurnal_var, params, q)
        real(dp), intent(in) :: y(:), daily_var(:), diurnal_var(:)
        type(mcsgarch_params_t), intent(in) :: params
        real(dp), intent(out) :: q(:)
        real(dp), allocatable :: e(:)
        integer :: t

        if (size(q) /= size(y)) error stop "mcsgarch_gjr_filter: q has wrong size"
        if (.not. params_valid(params, model_gjr)) error stop "mcsgarch_gjr_filter: invalid parameters"
        allocate(e(size(y)))
        e = y / sqrt(max(daily_var, min_var) * max(diurnal_var, min_var))
        q(1) = max(params%omega / max(1.0_dp - mcsgarch_persist(params, model_gjr), 1.0e-8_dp), variance(e))
        q(1) = max(q(1), min_var)
        do t = 2, size(y)
            q(t) = max(params%omega + params%alpha*e(t - 1)**2 + &
                       params%gamma*merge(e(t - 1)**2, 0.0_dp, e(t - 1) < 0.0_dp) + &
                       params%beta*q(t - 1), min_var)
        end do
        deallocate(e)
    end subroutine mcsgarch_gjr_filter

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
        if (.not. params_valid(params, model_sym)) error stop "mcsgarch_simulate: invalid parameters"
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
        real(dp) :: h, loss, nu, xi
        integer :: t

        call unpack_params(p, obj_model, obj_dist, params, nu, xi)
        if (.not. params_valid(params, obj_model)) then
            mcsgarch_nll_from_p = huge(1.0_dp) / 10.0_dp
            return
        end if
        allocate(q(size(obj_y)))
        select case (obj_model)
        case (model_nagarch)
            call mcsgarch_nagarch_filter(obj_y, obj_daily_var, obj_diurnal_var, params, q)
        case (model_gjr)
            call mcsgarch_gjr_filter(obj_y, obj_daily_var, obj_diurnal_var, params, q)
        case default
            call mcsgarch_filter(obj_y, obj_daily_var, obj_diurnal_var, params, q)
        end select
        loss = 0.0_dp
        do t = 1, size(obj_y)
            h = max(obj_daily_var(t) * obj_diurnal_var(t) * q(t), min_var)
            loss = loss + innovation_nll(obj_y(t), h, obj_dist, nu, xi)
        end do
        mcsgarch_nll_from_p = loss / real(size(obj_y), dp)
        deallocate(q)
    end function mcsgarch_nll_from_p

    ! Pack constrained model parameters into unconstrained optimizer variables.
    subroutine pack_params(params, model, dist, p)
        type(mcsgarch_params_t), intent(in) :: params
        integer, intent(in) :: model, dist
        real(dp), intent(out) :: p(:)
        real(dp) :: rest, a, b, ghalf, theta

        p(1) = log(max(params%omega, min_var))
        theta = params%theta
        if (model == model_nagarch) then
            a = max(params%alpha * (1.0_dp + theta**2), min_var)
        else
            a = max(params%alpha, min_var)
        end if
        b = max(params%beta, min_var)
        if (model == model_gjr) then
            ghalf = max(0.5_dp*params%gamma, min_var)
            rest = max(persist_max - a - b - ghalf, min_var)
        else
            ghalf = 0.0_dp
            rest = max(persist_max - a - b, min_var)
        end if
        p(2) = log(a / rest)
        p(3) = log(b / rest)
        if (model == model_nagarch) p(4) = theta
        if (model == model_gjr) p(4) = log(ghalf / rest)
        if (dist == dist_t) p(model_param_count(model) + 1) = log((8.0_dp - 2.0_dp) / (100.0_dp - 8.0_dp))
        if (dist == dist_fs_skewt) then
            p(model_param_count(model) + 1) = log((8.0_dp - 2.0_dp) / (100.0_dp - 8.0_dp))
            p(model_param_count(model) + 2) = 0.0_dp
        end if
    end subroutine pack_params

    ! Unpack optimizer variables into constrained model parameters.
    subroutine unpack_params(p, model, dist, params, nu, xi)
        real(dp), intent(in) :: p(:)
        integer, intent(in) :: model, dist
        type(mcsgarch_params_t), intent(out) :: params
        real(dp), intent(out) :: nu, xi
        real(dp) :: ea, eb, eg, den, theta, alpha_effective

        params%omega = exp(max(min(p(1), 50.0_dp), -50.0_dp))
        ea = exp(max(min(p(2), 50.0_dp), -50.0_dp))
        eb = exp(max(min(p(3), 50.0_dp), -50.0_dp))
        if (model == model_gjr) then
            eg = exp(max(min(p(4), 50.0_dp), -50.0_dp))
        else
            eg = 0.0_dp
        end if
        den = 1.0_dp + ea + eb + eg
        alpha_effective = persist_max * ea / den
        params%beta = persist_max * eb / den
        params%gamma = 0.0_dp
        if (model == model_nagarch) then
            theta = max(min(p(4), 20.0_dp), -20.0_dp)
            params%theta = theta
            params%alpha = alpha_effective / (1.0_dp + theta**2)
        else if (model == model_gjr) then
            params%theta = 0.0_dp
            params%alpha = alpha_effective
            params%gamma = 2.0_dp * persist_max * eg / den
        else
            params%theta = 0.0_dp
            params%alpha = alpha_effective
        end if
        if (dist == dist_t .or. dist == dist_fs_skewt) then
            nu = 2.0_dp + 98.0_dp / (1.0_dp + exp(-p(model_param_count(model) + 1)))
        else
            nu = 0.0_dp
        end if
        if (dist == dist_fs_skewt) then
            xi = exp(max(min(p(model_param_count(model) + 2), 20.0_dp), -20.0_dp))
        else
            xi = 0.0_dp
        end if
    end subroutine unpack_params

    ! Return the number of dynamic q_t parameters for a model variant.
    integer function model_param_count(model)
        integer, intent(in) :: model

        if (model == model_sym) then
            model_param_count = 3
        else
            model_param_count = 4
        end if
    end function model_param_count

    ! Negative log likelihood contribution for one zero-mean innovation.
    real(dp) function innovation_nll(y, h, dist, nu, xi)
        real(dp), intent(in) :: y, h, nu, xi
        integer, intent(in) :: dist
        real(dp) :: hh, z, pdf

        hh = max(h, min_var)
        if (dist == dist_t) then
            innovation_nll = log_gamma(0.5_dp*nu) - log_gamma(0.5_dp*(nu + 1.0_dp)) + &
                             0.5_dp*log(pi*(nu - 2.0_dp)) + 0.5_dp*log(hh) + &
                             0.5_dp*(nu + 1.0_dp)*log(1.0_dp + y**2 / ((nu - 2.0_dp)*hh))
        else if (dist == dist_fs_skewt) then
            z = y / sqrt(hh)
            pdf = max(pdf_fs_skewt(z, nu, xi), min_pdf)
            innovation_nll = 0.5_dp*log(hh) - log(pdf)
        else
            innovation_nll = log_sqrt_2pi + 0.5_dp*log(hh) + 0.5_dp*y**2/hh
        end if
    end function innovation_nll

    ! Return whether parameters satisfy positivity and persistence constraints.
    logical function params_valid(params, model)
        type(mcsgarch_params_t), intent(in) :: params
        integer, intent(in) :: model

        params_valid = params%omega > 0.0_dp .and. params%alpha >= 0.0_dp .and. &
                       params%beta >= 0.0_dp .and. params%gamma >= 0.0_dp .and. &
                       mcsgarch_persist(params, model) < 1.0_dp
    end function params_valid

    ! Return the persistence measure for the selected q_t dynamics.
    real(dp) function mcsgarch_persist(params, model)
        type(mcsgarch_params_t), intent(in) :: params
        integer, intent(in) :: model

        select case (model)
        case (model_nagarch)
            mcsgarch_persist = params%alpha * (1.0_dp + params%theta**2) + params%beta
        case (model_gjr)
            mcsgarch_persist = params%alpha + 0.5_dp*params%gamma + params%beta
        case default
            mcsgarch_persist = params%alpha + params%beta
        end select
    end function mcsgarch_persist

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

    ! Expected remaining-session variance from bar k given current GARCH state q_k.
    !
    ! Sums daily_var_d * diurnal_var_future(h) * E[q_{k+h} | F_k] over h = 1 to n_future,
    ! where E[q_{k+h} | F_k] = q_bar + persist^h * (q_k - q_bar)
    ! and q_bar = omega / (1 - persist).
    pure real(dp) function mcsgarch_remaining_var(params, persist, q_k, daily_var_d, diurnal_var_future)
        type(mcsgarch_params_t), intent(in) :: params
        real(dp), intent(in) :: persist, q_k, daily_var_d
        real(dp), intent(in) :: diurnal_var_future(:)
        integer :: h, n_future
        real(dp) :: q_bar, q_fcast, ppower

        n_future = size(diurnal_var_future)
        if (n_future == 0) then
            mcsgarch_remaining_var = 0.0_dp
            return
        end if
        q_bar = max(params%omega / max(1.0_dp - persist, 1.0e-10_dp), min_var)
        ppower = 1.0_dp
        mcsgarch_remaining_var = 0.0_dp
        do h = 1, n_future
            ppower = ppower * persist
            q_fcast = max(q_bar + ppower * (q_k - q_bar), min_var)
            mcsgarch_remaining_var = mcsgarch_remaining_var + diurnal_var_future(h) * q_fcast
        end do
        mcsgarch_remaining_var = daily_var_d * mcsgarch_remaining_var
    end function mcsgarch_remaining_var

end module garch_mcsgarch_mod

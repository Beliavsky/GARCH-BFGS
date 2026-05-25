module garch_split_fit_dist_mod
    use kind_mod, only: dp
    use math_const_mod, only: log_sqrt_2pi, sqrt2
    use garch_types_mod, only: garch_params_t
    use distributions_mod, only: pdf_normal, pdf_t, pdf_ged, pdf_logistic, pdf_laplace, pdf_sech, &
        pdf_nig, pdf_fs_skewt
    use garch_fit_dist_mod, only: dist_param_count
    use bfgs_mod, only: bfgs_minimize
    implicit none
    private

    public :: split_garch_dist_fit_result_t, fit_split_garch_dist_model, fit_split_garch_range_dist_model
    public :: split_model_param_count, split_garch_persist, split_garch_vol_forecast

    real(dp), parameter :: min_h = 1.0e-12_dp
    real(dp), parameter :: min_pdf = 1.0e-300_dp
    real(dp), parameter :: trading_days = 252.0_dp

    type :: split_garch_dist_fit_result_t
        character(len=16) :: model = ""
        character(len=10) :: dist = ""
        type(garch_params_t) :: params
        real(dp) :: co_frac = 0.0_dp
        real(dp) :: shape = 0.0_dp
        real(dp) :: shape2 = 0.0_dp
        real(dp) :: nll = huge(1.0_dp)
        real(dp) :: persist = 0.0_dp
        integer :: niter = 0
        logical :: converged = .false.
    end type split_garch_dist_fit_result_t

    real(dp), allocatable, save :: obj_cc(:), obj_co(:), obj_oc(:), obj_range(:)
    character(len=16), save :: obj_model = ""
    character(len=10), save :: obj_dist = ""
    integer, save :: obj_np = 0
    logical, save :: obj_use_range = .false.

contains

    subroutine fit_split_garch_dist_model(model_name, dist_name, cc, co, oc, max_iter, gtol, result)
        character(len=*), intent(in) :: model_name, dist_name
        real(dp), intent(in) :: cc(:), co(:), oc(:), gtol
        integer, intent(in) :: max_iter
        type(split_garch_dist_fit_result_t), intent(out) :: result
        real(dp), allocatable :: p(:), grad(:)
        real(dp) :: f_best
        integer :: niter
        logical :: converged

        if (size(cc) /= size(co) .or. size(cc) /= size(oc)) error stop "fit_split_garch_dist_model: inconsistent sizes"
        obj_use_range = .false.
        obj_model = canonical_split_model(model_name)
        obj_dist = canonical_dist(dist_name)
        obj_np = split_model_param_count(obj_model) + dist_param_count(obj_dist)
        if (allocated(obj_cc)) deallocate(obj_cc, obj_co, obj_oc)
        if (allocated(obj_range)) deallocate(obj_range)
        allocate(obj_cc(size(cc)), obj_co(size(co)), obj_oc(size(oc)))
        obj_cc = cc
        obj_co = co
        obj_oc = oc
        allocate(p(obj_np), grad(obj_np))

        call initial_params(obj_model, obj_dist, p)
        call bfgs_minimize(split_dist_obj, p, obj_np, max_iter, gtol, f_best, niter, converged)
        call unpack_all(obj_model, obj_dist, p, result%params, result%co_frac, result%shape, result%shape2)
        call split_dist_obj(p, obj_np, result%nll, grad)
        result%model = obj_model
        result%dist = obj_dist
        result%persist = split_garch_persist(obj_model, result%params)
        result%niter = niter
        result%converged = converged
        deallocate(p, grad)
    end subroutine fit_split_garch_dist_model

    subroutine fit_split_garch_range_dist_model(model_name, dist_name, cc, co, oc, range_var, max_iter, gtol, result)
        character(len=*), intent(in) :: model_name, dist_name
        real(dp), intent(in) :: cc(:), co(:), oc(:), range_var(:), gtol
        integer, intent(in) :: max_iter
        type(split_garch_dist_fit_result_t), intent(out) :: result
        real(dp), allocatable :: p(:), grad(:)
        real(dp) :: f_best
        integer :: niter
        logical :: converged

        if (size(cc) /= size(co) .or. size(cc) /= size(oc) .or. size(cc) /= size(range_var)) &
            error stop "fit_split_garch_range_dist_model: inconsistent sizes"
        obj_use_range = .true.
        obj_model = canonical_split_model(model_name)
        obj_dist = canonical_dist(dist_name)
        obj_np = split_model_param_count(obj_model) + dist_param_count(obj_dist)
        if (allocated(obj_cc)) deallocate(obj_cc, obj_co, obj_oc)
        if (allocated(obj_range)) deallocate(obj_range)
        allocate(obj_cc(size(cc)), obj_co(size(co)), obj_oc(size(oc)), obj_range(size(range_var)))
        obj_cc = cc
        obj_co = co
        obj_oc = oc
        obj_range = max(range_var, min_h)
        allocate(p(obj_np), grad(obj_np))

        call initial_params(obj_model, obj_dist, p)
        call bfgs_minimize(split_dist_obj, p, obj_np, max_iter, gtol, f_best, niter, converged)
        call unpack_all(obj_model, obj_dist, p, result%params, result%co_frac, result%shape, result%shape2)
        call split_dist_obj(p, obj_np, result%nll, grad)
        result%model = obj_model
        result%dist = obj_dist
        result%persist = split_garch_persist(obj_model, result%params)
        result%niter = niter
        result%converged = converged
        deallocate(p, grad)
    end subroutine fit_split_garch_range_dist_model

    subroutine split_dist_obj(p, np, f, g)
        integer, intent(in) :: np
        real(dp), intent(in) :: p(np)
        real(dp), intent(out) :: f, g(np)
        real(dp), allocatable :: pp(:), pm(:)
        real(dp) :: step, fp, fm
        integer :: i

        f = split_dist_value(p, np)
        allocate(pp(np), pm(np))
        do i = 1, np
            step = 1.0e-5_dp * max(1.0_dp, abs(p(i)))
            pp = p
            pm = p
            pp(i) = pp(i) + step
            pm(i) = pm(i) - step
            fp = split_dist_value(pp, np)
            fm = split_dist_value(pm, np)
            g(i) = (fp - fm) / (2.0_dp*step)
        end do
        deallocate(pp, pm)
    end subroutine split_dist_obj

    real(dp) function split_dist_value(p, np)
        integer, intent(in) :: np
        real(dp), intent(in) :: p(np)
        type(garch_params_t) :: params
        real(dp) :: co_frac, shape, shape2, persist, h, h_co, h_oc, zco, zoc
        integer :: t

        if (np /= obj_np) then
            split_dist_value = huge(1.0_dp)
            return
        end if
        call unpack_all(obj_model, obj_dist, p, params, co_frac, shape, shape2)
        persist = split_garch_persist(obj_model, params)
        h = initial_variance(obj_cc, params, persist)
        split_dist_value = 0.0_dp
        do t = 1, size(obj_cc)
            if (h <= 0.0_dp .or. h /= h .or. h > 1.0e100_dp) then
                split_dist_value = huge(1.0_dp)
                return
            end if
            h = max(h, min_h)
            h_co = max(co_frac*h, min_h)
            h_oc = max((1.0_dp - co_frac)*h, min_h)
            zco = obj_co(t) / sqrt(h_co)
            zoc = obj_oc(t) / sqrt(h_oc)
            split_dist_value = split_dist_value - log(max(innovation_pdf(obj_dist, zco, shape, shape2), min_pdf)) &
                + 0.5_dp*log(h_co) - log(max(innovation_pdf(obj_dist, zoc, shape, shape2), min_pdf)) &
                + 0.5_dp*log(h_oc)
            if (obj_use_range) split_dist_value = split_dist_value + log(h_oc) + obj_range(t)/h_oc
            h = next_variance(obj_cc(t), h, params, obj_model)
        end do
        split_dist_value = split_dist_value / real(size(obj_cc), dp)
        if (split_dist_value /= split_dist_value .or. split_dist_value > 1.0e29_dp) split_dist_value = 1.0e30_dp
    end function split_dist_value

    subroutine initial_params(model_name, dist_name, p)
        character(len=*), intent(in) :: model_name, dist_name
        real(dp), intent(out) :: p(:)
        type(garch_params_t) :: params
        real(dp) :: co_frac
        integer :: k

        params = garch_params_t()
        params%omega = max(sum(obj_cc**2) / real(size(obj_cc), dp), 1.0e-8_dp) * 0.03_dp
        params%alpha = 0.07_dp
        params%beta = 0.90_dp
        params%gamma = 0.05_dp
        params%theta = 0.5_dp
        params%twist = 0.0_dp
        co_frac = max(min(sum(obj_co**2) / max(sum(obj_co**2) + sum(obj_oc**2), min_h), 0.95_dp), 0.05_dp)
        call pack_split_params(model_name, params, co_frac, p)
        k = split_model_param_count(model_name) + 1
        call pack_dist_params(dist_name, p(k:))
    end subroutine initial_params

    subroutine unpack_all(model_name, dist_name, p, params, co_frac, shape, shape2)
        character(len=*), intent(in) :: model_name, dist_name
        real(dp), intent(in) :: p(:)
        type(garch_params_t), intent(out) :: params
        real(dp), intent(out) :: co_frac, shape, shape2
        integer :: k

        call unpack_split_params(model_name, p, params, co_frac)
        k = split_model_param_count(model_name) + 1
        call unpack_dist_params(dist_name, p(k:), shape, shape2)
    end subroutine unpack_all

    integer function split_model_param_count(model_name)
        character(len=*), intent(in) :: model_name
        select case (trim(model_name))
        case ("SPLIT_SYMM")
            split_model_param_count = 4
        case ("SPLIT_NAGARCH", "SPLIT_GJR")
            split_model_param_count = 5
        case ("SPLIT_FGTWIST")
            split_model_param_count = 6
        case default
            split_model_param_count = 0
        end select
    end function split_model_param_count

    subroutine pack_split_params(model_name, params, co_frac, p)
        character(len=*), intent(in) :: model_name
        type(garch_params_t), intent(in) :: params
        real(dp), intent(in) :: co_frac
        real(dp), intent(out) :: p(:)
        real(dp) :: slack, aa, bg, moment

        p(1) = log(max(params%omega, min_h))
        select case (trim(model_name))
        case ("SPLIT_SYMM")
            slack = max(1.0_dp - params%alpha - params%beta, 1.0e-8_dp)
            p(2) = log(max(params%alpha, min_h) / slack)
            p(3) = log(max(params%beta, min_h) / slack)
            p(4) = logit(co_frac)
        case ("SPLIT_NAGARCH")
            moment = 1.0_dp + params%theta**2
            aa = params%alpha * moment
            slack = max(1.0_dp - aa - params%beta, 1.0e-8_dp)
            p(2) = log(max(aa, min_h) / slack)
            p(3) = log(max(params%beta, min_h) / slack)
            p(4) = params%theta
            p(5) = logit(co_frac)
        case ("SPLIT_GJR")
            bg = 0.5_dp * params%gamma
            slack = max(1.0_dp - params%alpha - bg - params%beta, 1.0e-8_dp)
            p(2) = log(max(params%alpha, min_h) / slack)
            p(3) = log(max(bg, min_h) / slack)
            p(4) = log(max(params%beta, min_h) / slack)
            p(5) = logit(co_frac)
        case ("SPLIT_FGTWIST")
            moment = fgarch_twist_moment_local(params%theta, params%twist)
            aa = params%alpha * moment
            slack = max(1.0_dp - aa - params%beta, 1.0e-8_dp)
            p(2) = log(max(aa, min_h) / slack)
            p(3) = log(max(params%beta, min_h) / slack)
            p(4) = params%theta
            p(5) = params%twist
            p(6) = logit(co_frac)
        end select
    end subroutine pack_split_params

    subroutine unpack_split_params(model_name, p, params, co_frac)
        character(len=*), intent(in) :: model_name
        real(dp), intent(in) :: p(:)
        type(garch_params_t), intent(out) :: params
        real(dp), intent(out) :: co_frac
        real(dp) :: e2, e3, e4, s, aa, moment

        params = garch_params_t()
        params%omega = exp(clamp(p(1), -60.0_dp, 20.0_dp))
        select case (trim(model_name))
        case ("SPLIT_SYMM")
            call simplex2(p(2), p(3), params%alpha, params%beta)
            co_frac = logistic(p(4))
        case ("SPLIT_NAGARCH")
            params%theta = p(4)
            moment = 1.0_dp + params%theta**2
            call simplex2(p(2), p(3), aa, params%beta)
            params%alpha = aa / moment
            co_frac = logistic(p(5))
        case ("SPLIT_GJR")
            e2 = exp(clamp(p(2), -40.0_dp, 40.0_dp))
            e3 = exp(clamp(p(3), -40.0_dp, 40.0_dp))
            e4 = exp(clamp(p(4), -40.0_dp, 40.0_dp))
            s = 1.0_dp + e2 + e3 + e4
            params%alpha = e2 / s
            params%gamma = 2.0_dp * e3 / s
            params%beta = e4 / s
            co_frac = logistic(p(5))
        case ("SPLIT_FGTWIST")
            params%theta = p(4)
            params%twist = p(5)
            moment = fgarch_twist_moment_local(params%theta, params%twist)
            call simplex2(p(2), p(3), aa, params%beta)
            params%alpha = aa / moment
            co_frac = logistic(p(6))
        end select
    end subroutine unpack_split_params

    subroutine pack_dist_params(dist_name, p)
        character(len=*), intent(in) :: dist_name
        real(dp), intent(out) :: p(:)
        if (size(p) == 0) return
        select case (trim(dist_name))
        case ("T")
            p(1) = pack_bounded(8.0_dp, 2.01_dp, 100.0_dp)
        case ("GED")
            p(1) = log(1.5_dp)
        case ("NIG")
            p(1) = pack_bounded(3.0_dp, 0.1_dp, 20.0_dp)
        case ("FS_SKEWT")
            p(1) = pack_bounded(8.0_dp, 2.01_dp, 100.0_dp)
            p(2) = 0.0_dp
        case default
            p = 0.0_dp
        end select
    end subroutine pack_dist_params

    subroutine unpack_dist_params(dist_name, p, shape, shape2)
        character(len=*), intent(in) :: dist_name
        real(dp), intent(in) :: p(:)
        real(dp), intent(out) :: shape, shape2
        shape = 0.0_dp
        shape2 = 0.0_dp
        select case (trim(dist_name))
        case ("T")
            shape = unpack_bounded(p(1), 2.01_dp, 100.0_dp)
        case ("GED")
            shape = exp(clamp(p(1), -40.0_dp, 40.0_dp))
        case ("NIG")
            shape = unpack_bounded(p(1), 0.1_dp, 20.0_dp)
        case ("FS_SKEWT")
            shape = unpack_bounded(p(1), 2.01_dp, 100.0_dp)
            shape2 = exp(clamp(p(2), -40.0_dp, 40.0_dp))
        end select
    end subroutine unpack_dist_params

    real(dp) function innovation_pdf(dist_name, z, shape, shape2)
        character(len=*), intent(in) :: dist_name
        real(dp), intent(in) :: z, shape, shape2
        select case (trim(dist_name))
        case ("NORMAL")
            innovation_pdf = pdf_normal(z)
        case ("T")
            innovation_pdf = pdf_t(z, shape)
        case ("SECH")
            innovation_pdf = pdf_sech(z)
        case ("GED")
            innovation_pdf = pdf_ged(z, shape)
        case ("LAPLACE")
            innovation_pdf = pdf_laplace(z)
        case ("LOGISTIC")
            innovation_pdf = pdf_logistic(z)
        case ("NIG")
            innovation_pdf = pdf_nig(z, shape)
        case ("FS_SKEWT")
            innovation_pdf = pdf_fs_skewt(z, shape, shape2)
        case default
            innovation_pdf = min_pdf
        end select
    end function innovation_pdf

    real(dp) function split_garch_persist(model_name, params)
        character(len=*), intent(in) :: model_name
        type(garch_params_t), intent(in) :: params
        select case (trim(model_name))
        case ("SPLIT_SYMM")
            split_garch_persist = params%alpha + params%beta
        case ("SPLIT_NAGARCH")
            split_garch_persist = params%alpha*(1.0_dp + params%theta**2) + params%beta
        case ("SPLIT_GJR")
            split_garch_persist = params%alpha + 0.5_dp*params%gamma + params%beta
        case ("SPLIT_FGTWIST")
            split_garch_persist = params%alpha*fgarch_twist_moment_local(params%theta, params%twist) + params%beta
        case default
            split_garch_persist = 0.0_dp
        end select
    end function split_garch_persist

    subroutine split_garch_vol_forecast(cc, model_name, params, persist, vol)
        real(dp), intent(in) :: cc(:), persist
        character(len=*), intent(in) :: model_name
        type(garch_params_t), intent(in) :: params
        real(dp), intent(out) :: vol(:)
        real(dp) :: h
        integer :: t

        h = initial_variance(cc, params, persist)
        do t = 1, size(cc)
            h = max(h, min_h)
            vol(t) = sqrt(trading_days*h) * 100.0_dp
            h = next_variance(cc(t), h, params, model_name)
        end do
    end subroutine split_garch_vol_forecast

    real(dp) function initial_variance(y, params, persist)
        real(dp), intent(in) :: y(:), persist
        type(garch_params_t), intent(in) :: params
        initial_variance = max(params%omega / max(1.0_dp - persist, 1.0e-8_dp), &
            sum(y**2) / real(size(y), dp), min_h)
    end function initial_variance

    real(dp) function next_variance(y, h, params, model_name)
        real(dp), intent(in) :: y, h
        type(garch_params_t), intent(in) :: params
        character(len=*), intent(in) :: model_name
        real(dp) :: sqrth, r, ind, z, q

        sqrth = sqrt(max(h, min_h))
        select case (trim(model_name))
        case ("SPLIT_SYMM")
            next_variance = params%omega + params%alpha*y**2 + params%beta*h
        case ("SPLIT_NAGARCH")
            r = y - params%theta*sqrth
            next_variance = params%omega + params%alpha*r**2 + params%beta*h
        case ("SPLIT_GJR")
            ind = merge(1.0_dp, 0.0_dp, y < 0.0_dp)
            next_variance = params%omega + (params%alpha + params%gamma*ind)*y**2 + params%beta*h
        case ("SPLIT_FGTWIST")
            z = y / sqrth
            q = abs(z - params%theta) - params%twist*(z - params%theta)
            next_variance = params%omega + params%alpha*h*q**2 + params%beta*h
        case default
            next_variance = h
        end select
        next_variance = max(next_variance, min_h)
    end function next_variance

    subroutine simplex2(pa, pb, a, b)
        real(dp), intent(in) :: pa, pb
        real(dp), intent(out) :: a, b
        real(dp) :: ea, eb, s
        ea = exp(clamp(pa, -40.0_dp, 40.0_dp))
        eb = exp(clamp(pb, -40.0_dp, 40.0_dp))
        s = 1.0_dp + ea + eb
        a = ea / s
        b = eb / s
    end subroutine simplex2

    real(dp) function logistic(x)
        real(dp), intent(in) :: x
        logistic = 1.0_dp / (1.0_dp + exp(-clamp(x, -40.0_dp, 40.0_dp)))
        logistic = min(max(logistic, 1.0e-5_dp), 1.0_dp - 1.0e-5_dp)
    end function logistic

    real(dp) function logit(x)
        real(dp), intent(in) :: x
        real(dp) :: xx
        xx = min(max(x, 1.0e-5_dp), 1.0_dp - 1.0e-5_dp)
        logit = log(xx / (1.0_dp - xx))
    end function logit

    real(dp) function pack_bounded(x, lo, hi)
        real(dp), intent(in) :: x, lo, hi
        real(dp) :: xx
        xx = min(max(x, lo + 1.0e-8_dp), hi - 1.0e-8_dp)
        pack_bounded = log((xx - lo) / (hi - xx))
    end function pack_bounded

    real(dp) function unpack_bounded(q, lo, hi)
        real(dp), intent(in) :: q, lo, hi
        unpack_bounded = lo + (hi - lo) / (1.0_dp + exp(-clamp(q, -40.0_dp, 40.0_dp)))
    end function unpack_bounded

    real(dp) function clamp(x, lo, hi)
        real(dp), intent(in) :: x, lo, hi
        clamp = min(max(x, lo), hi)
    end function clamp

    real(dp) function fgarch_twist_moment_local(theta, twist)
        real(dp), intent(in) :: theta, twist
        real(dp) :: Phi, ph, A, B, one_minus, one_plus
        Phi = 0.5_dp * (1.0_dp + erf(theta / sqrt2))
        ph  = exp(-0.5_dp * theta**2 - log_sqrt_2pi)
        A = (1.0_dp + theta**2) * (1.0_dp - Phi) - theta * ph
        B = (1.0_dp + theta**2) * Phi + theta * ph
        one_minus = 1.0_dp - twist
        one_plus  = 1.0_dp + twist
        fgarch_twist_moment_local = max(one_minus**2 * A + one_plus**2 * B, 1.0e-8_dp)
    end function fgarch_twist_moment_local

    character(len=16) function canonical_split_model(model_name)
        character(len=*), intent(in) :: model_name
        canonical_split_model = adjustl(model_name)
    end function canonical_split_model

    character(len=10) function canonical_dist(dist_name)
        character(len=*), intent(in) :: dist_name
        canonical_dist = adjustl(dist_name)
    end function canonical_dist

end module garch_split_fit_dist_mod

module garch_fit_dist_mod
    use kind_mod, only: dp
    use garch_types_mod, only: garch_params_t
    use distributions_mod, only: pdf_normal, pdf_t, pdf_ged, pdf_logistic, pdf_laplace, pdf_sech, pdf_nig, pdf_fs_skewt
    use garch_fit_mod, only: fit_symm_garch, fit_nagarch, fit_gjr, fit_gjr_signed, fit_egarch, fit_qgarch
    use stats_mod, only: variance
    use bfgs_mod, only: bfgs_minimize
    implicit none
    private

    integer, parameter :: dist_normal = 1
    integer, parameter :: dist_t = 2
    integer, parameter :: dist_sech = 3
    integer, parameter :: dist_ged = 4
    integer, parameter :: dist_laplace = 5
    integer, parameter :: dist_logistic = 6
    integer, parameter :: dist_nig = 7
    integer, parameter :: dist_fs_skewt = 8
    real(dp), parameter :: min_h = 1.0e-12_dp
    real(dp), parameter :: min_pdf = 1.0e-300_dp
    real(dp), parameter :: persist_max = 0.999_dp

    type, public :: garch_dist_fit_result_t
        character(len=16) :: model = ""
        character(len=10) :: dist = ""
        type(garch_params_t) :: params
        real(dp) :: shape = 0.0_dp
        real(dp) :: shape2 = 0.0_dp
        real(dp) :: nll = huge(1.0_dp)
        real(dp) :: persist = 0.0_dp
        integer :: niter = 0
        logical :: converged = .false.
    end type garch_dist_fit_result_t

    real(dp), allocatable, save :: obj_y(:)
    character(len=16), save :: obj_model = ""
    integer, save :: obj_dist = dist_normal
    integer, save :: obj_np = 0

    public :: fit_garch_dist_model, garch_dist_oos_nll, garch_dist_persist, garch_dist_variance_path
    public :: model_param_count, dist_param_count, dist_nu_value, dist_xi_value, dist_alpha_value

contains

    subroutine fit_garch_dist_model(model_name, dist_name, y, max_iter, gtol, result, &
                                    start_params, start_shape, start_shape2)
        character(len=*), intent(in) :: model_name, dist_name
        real(dp), intent(in) :: y(:), gtol
        integer, intent(in) :: max_iter
        type(garch_dist_fit_result_t), intent(out) :: result
        type(garch_params_t), intent(in), optional :: start_params
        real(dp), intent(in), optional :: start_shape, start_shape2
        type(garch_params_t) :: params0
        real(dp) :: shape0, shape20
        real(dp), allocatable :: p(:), grad(:)
        real(dp) :: f_best
        integer :: niter
        logical :: converged

        if (present(start_params)) then
            params0 = start_params
        else
            call initial_params(model_name, y, max_iter, gtol, params0)
        end if
        shape0 = default_shape(dist_id_from_name(dist_name))
        shape20 = default_shape2(dist_id_from_name(dist_name))
        if (present(start_shape)) shape0 = start_shape
        if (present(start_shape2)) shape20 = start_shape2
        obj_model = canonical_model(model_name)
        obj_dist = dist_id_from_name(dist_name)
        obj_np = model_npar(obj_model) + dist_nshape(obj_dist)
        if (allocated(obj_y)) deallocate(obj_y)
        allocate(obj_y(size(y)))
        obj_y = y
        allocate(p(obj_np), grad(obj_np))

        call pack_params(obj_model, obj_dist, params0, shape0, shape20, p)
        call bfgs_minimize(dist_obj, p, obj_np, max_iter, gtol, f_best, niter, converged)
        call unpack_params(obj_model, obj_dist, p, result%params, result%shape, result%shape2)
        call dist_obj(p, obj_np, result%nll, grad)

        result%model = obj_model
        result%dist = canonical_dist(dist_name)
        result%persist = garch_dist_persist(obj_model, result%params)
        result%niter = niter
        result%converged = converged

        deallocate(p, grad)
    end subroutine fit_garch_dist_model

    subroutine dist_obj(p, np, f, g)
        integer, intent(in) :: np
        real(dp), intent(in) :: p(np)
        real(dp), intent(out) :: f, g(np)
        real(dp) :: step, fp, fm
        real(dp), allocatable :: pp(:), pm(:)
        integer :: j

        f = nll_from_p(p, np)
        allocate(pp(np), pm(np))
        do j = 1, np
            step = 1.0e-5_dp * max(1.0_dp, abs(p(j)))
            pp = p
            pm = p
            pp(j) = pp(j) + step
            pm(j) = pm(j) - step
            fp = nll_from_p(pp, np)
            fm = nll_from_p(pm, np)
            g(j) = (fp - fm) / (2.0_dp * step)
        end do
        deallocate(pp, pm)
    end subroutine dist_obj

    real(dp) function nll_from_p(p, np)
        integer, intent(in) :: np
        real(dp), intent(in) :: p(np)
        type(garch_params_t) :: params
        real(dp) :: shape, shape2

        if (np /= obj_np) then
            nll_from_p = huge(1.0_dp)
            return
        end if
        call unpack_params(obj_model, obj_dist, p, params, shape, shape2)
        nll_from_p = nll_model_dist(obj_model, obj_dist, obj_y, params, shape, shape2)
    end function nll_from_p

    real(dp) function garch_dist_oos_nll(model_name, dist_name, y_full, ntrain, ntest, result)
        character(len=*), intent(in) :: model_name, dist_name
        real(dp), intent(in) :: y_full(:)
        integer, intent(in) :: ntrain, ntest
        type(garch_dist_fit_result_t), intent(in) :: result
        real(dp), allocatable :: h(:)
        real(dp) :: z, dens, loss
        integer :: t, dist_id

        dist_id = dist_id_from_name(dist_name)
        allocate(h(size(y_full)))
        call variance_path(canonical_model(model_name), y_full, result%params, h)
        loss = 0.0_dp
        do t = ntrain + 1, ntrain + ntest
            z = y_full(t) / sqrt(max(h(t), min_h))
            dens = max(innovation_pdf(dist_id, z, result%shape, result%shape2), min_pdf)
            loss = loss - log(dens) + 0.5_dp*log(max(h(t), min_h))
        end do
        garch_dist_oos_nll = loss / real(ntest, dp)
        deallocate(h)
    end function garch_dist_oos_nll

    subroutine garch_dist_variance_path(model_name, y, params, h)
        character(len=*), intent(in) :: model_name
        real(dp), intent(in) :: y(:)
        type(garch_params_t), intent(in) :: params
        real(dp), intent(out) :: h(:)

        call variance_path(canonical_model(model_name), y, params, h)
    end subroutine garch_dist_variance_path

    real(dp) function nll_model_dist(model_name, dist_id, y, params, shape, shape2)
        character(len=*), intent(in) :: model_name
        integer, intent(in) :: dist_id
        real(dp), intent(in) :: y(:), shape, shape2
        type(garch_params_t), intent(in) :: params
        real(dp), allocatable :: h(:)
        real(dp) :: z, dens, loss
        integer :: t

        if (.not. params_valid(model_name, params)) then
            nll_model_dist = huge(1.0_dp) / 10.0_dp
            return
        end if
        allocate(h(size(y)))
        call variance_path(model_name, y, params, h)
        loss = 0.0_dp
        do t = 1, size(y)
            z = y(t) / sqrt(max(h(t), min_h))
            dens = max(innovation_pdf(dist_id, z, shape, shape2), min_pdf)
            loss = loss - log(dens) + 0.5_dp*log(max(h(t), min_h))
        end do
        nll_model_dist = loss / real(size(y), dp)
        deallocate(h)
    end function nll_model_dist

    subroutine variance_path(model_name, y, params, h)
        character(len=*), intent(in) :: model_name
        real(dp), intent(in) :: y(:)
        type(garch_params_t), intent(in) :: params
        real(dp), intent(out) :: h(:)
        real(dp) :: lh, z, c_eg, sqrth, ind, var0
        integer :: t

        var0 = max(variance(y), min_h)
        select case (model_name)
        case ("SYMM_GARCH", "GARCH")
            h(1) = max(params%omega / max(1.0_dp - params%alpha - params%beta, 1.0e-8_dp), var0)
            do t = 2, size(y)
                h(t) = max(params%omega + params%alpha*y(t - 1)**2 + params%beta*h(t - 1), min_h)
            end do
        case ("NAGARCH")
            h(1) = max(params%omega / max(1.0_dp - garch_dist_persist(model_name, params), 1.0e-8_dp), var0)
            do t = 2, size(y)
                sqrth = sqrt(max(h(t - 1), min_h))
                h(t) = max(params%omega + params%alpha*(y(t - 1) - params%theta*sqrth)**2 + &
                           params%beta*h(t - 1), min_h)
            end do
        case ("GJR_GARCH", "GJR")
            h(1) = max(params%omega / max(1.0_dp - garch_dist_persist(model_name, params), 1.0e-8_dp), var0)
            do t = 2, size(y)
                ind = merge(1.0_dp, 0.0_dp, y(t - 1) < 0.0_dp)
                h(t) = max(params%omega + (params%alpha + params%gamma*ind)*y(t - 1)**2 + &
                           params%beta*h(t - 1), min_h)
            end do
        case ("GJR_SIGNED")
            h(1) = max(params%omega / max(1.0_dp - garch_dist_persist(model_name, params), 1.0e-8_dp), var0)
            do t = 2, size(y)
                ind = sign(1.0_dp, y(t - 1))
                h(t) = max(params%omega + max(params%alpha + params%gamma*ind, 0.0_dp)*y(t - 1)**2 + &
                           params%beta*h(t - 1), min_h)
            end do
        case ("EGARCH")
            c_eg = sqrt(2.0_dp / acos(-1.0_dp))
            lh = params%omega / max(1.0_dp - params%beta, 1.0e-8_dp)
            do t = 1, size(y)
                h(t) = max(exp(lh), min_h)
                z = y(t) / sqrt(h(t))
                lh = params%omega + params%beta*lh + params%alpha*(abs(z) - c_eg) + params%gamma*z
            end do
        case ("QGARCH")
            h(1) = max((params%omega + params%alpha*params%theta**2) / &
                       max(1.0_dp - params%alpha - params%beta, 1.0e-8_dp), var0)
            do t = 2, size(y)
                h(t) = max(params%omega + params%alpha*(y(t - 1) - params%theta)**2 + &
                           params%beta*h(t - 1), min_h)
            end do
        case default
            error stop "variance_path: unsupported model"
        end select
        h = max(h, min_h)
    end subroutine variance_path

    subroutine initial_params(model_name, y, max_iter, gtol, params)
        character(len=*), intent(in) :: model_name
        real(dp), intent(in) :: y(:), gtol
        integer, intent(in) :: max_iter
        type(garch_params_t), intent(out) :: params
        real(dp) :: f_best
        integer :: niter
        logical :: converged

        select case (canonical_model(model_name))
        case ("SYMM_GARCH", "GARCH")
            call fit_symm_garch(y, max_iter, gtol, f_best, params, niter, converged)
        case ("NAGARCH")
            call fit_nagarch(y, max_iter, gtol, f_best, params, niter, converged)
        case ("GJR_GARCH", "GJR")
            call fit_gjr(y, max_iter, gtol, f_best, params, niter, converged)
        case ("GJR_SIGNED")
            call fit_gjr_signed(y, max_iter, gtol, f_best, params, niter, converged)
        case ("EGARCH")
            call fit_egarch(y, max_iter, gtol, f_best, params, niter, converged)
        case ("QGARCH")
            call fit_qgarch(y, max_iter, gtol, f_best, params, niter, converged)
        case default
            error stop "initial_params: unsupported model"
        end select
    end subroutine initial_params

    subroutine pack_params(model_name, dist_id, params, shape, shape2, p)
        character(len=*), intent(in) :: model_name
        integer, intent(in) :: dist_id
        type(garch_params_t), intent(in) :: params
        real(dp), intent(in) :: shape, shape2
        real(dp), intent(out) :: p(:)
        integer :: k

        select case (model_name)
        case ("SYMM_GARCH", "GARCH", "QGARCH")
            p(1) = log(max(params%omega, min_h))
            call pack_two_simplex(params%alpha, params%beta, p(2), p(3))
            k = 4
            if (model_name == "QGARCH") then
                p(4) = params%theta
                k = 5
            end if
        case ("NAGARCH")
            p(1) = log(max(params%omega, min_h))
            call pack_two_simplex(params%alpha*(1.0_dp + params%theta**2), params%beta, p(2), p(3))
            p(4) = params%theta
            k = 5
        case ("GJR_GARCH", "GJR")
            p(1) = log(max(params%omega, min_h))
            call pack_three_simplex(params%alpha, 0.5_dp*params%gamma, params%beta, p(2), p(3), p(4))
            k = 5
        case ("GJR_SIGNED")
            p(1) = log(max(params%omega, min_h))
            call pack_two_simplex(params%alpha, params%beta, p(2), p(3))
            p(4) = atanh_clip(params%gamma / 0.5_dp)
            k = 5
        case ("EGARCH")
            p(1) = params%omega
            p(2) = params%alpha
            p(3) = params%gamma
            p(4) = logit_clip(params%beta / persist_max)
            k = 5
        case default
            error stop "pack_params: unsupported model"
        end select
        if (dist_nshape(dist_id) >= 1) p(k) = pack_shape(dist_id, shape)
        if (dist_nshape(dist_id) >= 2) p(k + 1) = pack_shape2(dist_id, shape2)
    end subroutine pack_params

    subroutine unpack_params(model_name, dist_id, p, params, shape, shape2)
        character(len=*), intent(in) :: model_name
        integer, intent(in) :: dist_id
        real(dp), intent(in) :: p(:)
        type(garch_params_t), intent(out) :: params
        real(dp), intent(out) :: shape, shape2
        real(dp) :: comp1, comp2, comp3
        integer :: k

        params = garch_params_t()
        select case (model_name)
        case ("SYMM_GARCH", "GARCH")
            params%omega = exp(clamp(p(1), -60.0_dp, 20.0_dp))
            call unpack_two_simplex(p(2), p(3), params%alpha, params%beta)
            k = 4
        case ("QGARCH")
            params%omega = exp(clamp(p(1), -60.0_dp, 20.0_dp))
            call unpack_two_simplex(p(2), p(3), params%alpha, params%beta)
            params%theta = p(4)
            k = 5
        case ("NAGARCH")
            params%omega = exp(clamp(p(1), -60.0_dp, 20.0_dp))
            call unpack_two_simplex(p(2), p(3), comp1, params%beta)
            params%theta = p(4)
            params%alpha = comp1 / (1.0_dp + params%theta**2)
            k = 5
        case ("GJR_GARCH", "GJR")
            params%omega = exp(clamp(p(1), -60.0_dp, 20.0_dp))
            call unpack_three_simplex(p(2), p(3), p(4), comp1, comp2, comp3)
            params%alpha = comp1
            params%gamma = 2.0_dp*comp2
            params%beta = comp3
            k = 5
        case ("GJR_SIGNED")
            params%omega = exp(clamp(p(1), -60.0_dp, 20.0_dp))
            call unpack_two_simplex(p(2), p(3), params%alpha, params%beta)
            params%gamma = 0.5_dp*tanh(p(4))
            k = 5
        case ("EGARCH")
            params%omega = p(1)
            params%alpha = p(2)
            params%gamma = p(3)
            params%beta = persist_max*sigmoid(p(4))
            k = 5
        case default
            error stop "unpack_params: unsupported model"
        end select
        if (dist_nshape(dist_id) >= 1) then
            shape = unpack_shape(dist_id, p(k))
        else
            shape = 0.0_dp
        end if
        if (dist_nshape(dist_id) >= 2) then
            shape2 = unpack_shape2(dist_id, p(k + 1))
        else
            shape2 = 0.0_dp
        end if
    end subroutine unpack_params

    subroutine pack_two_simplex(alpha, beta, pa, pb)
        real(dp), intent(in) :: alpha, beta
        real(dp), intent(out) :: pa, pb
        real(dp) :: base
        base = max(persist_max - alpha - beta, 1.0e-6_dp)
        pa = log(max(alpha, 1.0e-8_dp) / base)
        pb = log(max(beta, 1.0e-8_dp) / base)
    end subroutine pack_two_simplex

    subroutine unpack_two_simplex(pa, pb, alpha, beta)
        real(dp), intent(in) :: pa, pb
        real(dp), intent(out) :: alpha, beta
        real(dp) :: ea, eb, denom
        ea = exp(clamp(pa, -40.0_dp, 40.0_dp))
        eb = exp(clamp(pb, -40.0_dp, 40.0_dp))
        denom = 1.0_dp + ea + eb
        alpha = persist_max * ea / denom
        beta = persist_max * eb / denom
    end subroutine unpack_two_simplex

    subroutine pack_three_simplex(a, ghalf, b, pa, pg, pb)
        real(dp), intent(in) :: a, ghalf, b
        real(dp), intent(out) :: pa, pg, pb
        real(dp) :: base
        base = max(persist_max - a - ghalf - b, 1.0e-6_dp)
        pa = log(max(a, 1.0e-8_dp) / base)
        pg = log(max(ghalf, 1.0e-8_dp) / base)
        pb = log(max(b, 1.0e-8_dp) / base)
    end subroutine pack_three_simplex

    subroutine unpack_three_simplex(pa, pg, pb, a, ghalf, b)
        real(dp), intent(in) :: pa, pg, pb
        real(dp), intent(out) :: a, ghalf, b
        real(dp) :: ea, eg, eb, denom
        ea = exp(clamp(pa, -40.0_dp, 40.0_dp))
        eg = exp(clamp(pg, -40.0_dp, 40.0_dp))
        eb = exp(clamp(pb, -40.0_dp, 40.0_dp))
        denom = 1.0_dp + ea + eg + eb
        a = persist_max * ea / denom
        ghalf = persist_max * eg / denom
        b = persist_max * eb / denom
    end subroutine unpack_three_simplex

    integer function model_npar(model_name)
        character(len=*), intent(in) :: model_name
        select case (model_name)
        case ("SYMM_GARCH", "GARCH")
            model_npar = 3
        case ("NAGARCH", "GJR_GARCH", "GJR", "GJR_SIGNED", "EGARCH", "QGARCH")
            model_npar = 4
        case default
            error stop "model_npar: unsupported model"
        end select
    end function model_npar

    integer function model_param_count(model_name)
        character(len=*), intent(in) :: model_name

        select case (trim(model_name))
        case ("SYMM_GARCH", "GARCH")
            model_param_count = 3
        case ("NAGARCH", "GJR_GARCH", "GJR", "GJR_SIGNED", "EGARCH", "QGARCH")
            model_param_count = 4
        case default
            model_param_count = 0
        end select
    end function model_param_count

    integer function dist_nshape(dist_id)
        integer, intent(in) :: dist_id
        select case (dist_id)
        case (dist_t, dist_ged, dist_nig)
            dist_nshape = 1
        case (dist_fs_skewt)
            dist_nshape = 2
        case default
            dist_nshape = 0
        end select
    end function dist_nshape

    integer function dist_param_count(dist_name)
        character(len=*), intent(in) :: dist_name

        select case (trim(dist_name))
        case ("T", "GED", "NIG")
            dist_param_count = 1
        case ("FS_SKEWT")
            dist_param_count = 2
        case default
            dist_param_count = 0
        end select
    end function dist_param_count

    character(len=8) function dist_nu_value(dist_name, shape)
        character(len=*), intent(in) :: dist_name
        real(dp), intent(in) :: shape

        select case (trim(dist_name))
        case ("FS_SKEWT")
            write(dist_nu_value, '(F8.3)') shape
        case ("T", "GED")
            write(dist_nu_value, '(F8.3)') shape
        case default
            dist_nu_value = "-"
        end select
    end function dist_nu_value

    character(len=8) function dist_xi_value(dist_name, shape2)
        character(len=*), intent(in) :: dist_name
        real(dp), intent(in) :: shape2

        select case (trim(dist_name))
        case ("FS_SKEWT")
            write(dist_xi_value, '(F8.3)') shape2
        case default
            dist_xi_value = "-"
        end select
    end function dist_xi_value

    character(len=10) function dist_alpha_value(dist_name, shape)
        character(len=*), intent(in) :: dist_name
        real(dp), intent(in) :: shape

        select case (trim(dist_name))
        case ("NIG")
            write(dist_alpha_value, '(F10.3)') shape
        case default
            dist_alpha_value = "-"
        end select
    end function dist_alpha_value

    real(dp) function default_shape(dist_id)
        integer, intent(in) :: dist_id
        select case (dist_id)
        case (dist_t)
            default_shape = 8.0_dp
        case (dist_ged)
            default_shape = 1.5_dp
        case (dist_nig)
            default_shape = 3.0_dp
        case (dist_fs_skewt)
            default_shape = 8.0_dp
        case default
            default_shape = 0.0_dp
        end select
    end function default_shape

    real(dp) function pack_shape(dist_id, shape)
        integer, intent(in) :: dist_id
        real(dp), intent(in) :: shape
        select case (dist_id)
        case (dist_t)
            pack_shape = log((max(shape, 2.01_dp) - 2.0_dp) / max(100.0_dp - min(shape, 99.9_dp), 1.0e-6_dp))
        case (dist_ged)
            pack_shape = log(max(shape, 1.0e-4_dp))
        case (dist_nig)
            pack_shape = log((max(shape, 0.11_dp) - 0.1_dp) / max(20.0_dp - min(shape, 19.9_dp), 1.0e-6_dp))
        case (dist_fs_skewt)
            pack_shape = log((max(shape, 2.01_dp) - 2.0_dp) / max(100.0_dp - min(shape, 99.9_dp), 1.0e-6_dp))
        case default
            pack_shape = 0.0_dp
        end select
    end function pack_shape

    real(dp) function unpack_shape(dist_id, p)
        integer, intent(in) :: dist_id
        real(dp), intent(in) :: p
        select case (dist_id)
        case (dist_t)
            unpack_shape = 2.0_dp + 98.0_dp*sigmoid(p)
        case (dist_ged)
            unpack_shape = exp(clamp(p, -8.0_dp, 4.0_dp))
        case (dist_nig)
            unpack_shape = 0.1_dp + 19.9_dp*sigmoid(p)
        case (dist_fs_skewt)
            unpack_shape = 2.0_dp + 98.0_dp*sigmoid(p)
        case default
            unpack_shape = 0.0_dp
        end select
    end function unpack_shape

    real(dp) function default_shape2(dist_id)
        integer, intent(in) :: dist_id
        select case (dist_id)
        case (dist_fs_skewt)
            default_shape2 = 1.0_dp
        case default
            default_shape2 = 0.0_dp
        end select
    end function default_shape2

    real(dp) function pack_shape2(dist_id, shape2)
        integer, intent(in) :: dist_id
        real(dp), intent(in) :: shape2
        select case (dist_id)
        case (dist_fs_skewt)
            pack_shape2 = log(max(shape2, 1.0e-6_dp))
        case default
            pack_shape2 = 0.0_dp
        end select
    end function pack_shape2

    real(dp) function unpack_shape2(dist_id, p)
        integer, intent(in) :: dist_id
        real(dp), intent(in) :: p
        select case (dist_id)
        case (dist_fs_skewt)
            unpack_shape2 = exp(clamp(p, -4.0_dp, 4.0_dp))
        case default
            unpack_shape2 = 0.0_dp
        end select
    end function unpack_shape2

    logical function params_valid(model_name, params)
        character(len=*), intent(in) :: model_name
        type(garch_params_t), intent(in) :: params
        params_valid = garch_dist_persist(model_name, params) < persist_max .and. &
                       garch_dist_persist(model_name, params) > -0.999_dp
        if (model_name /= "EGARCH") params_valid = params_valid .and. params%omega > 0.0_dp
    end function params_valid

    real(dp) function garch_dist_persist(model_name, params)
        character(len=*), intent(in) :: model_name
        type(garch_params_t), intent(in) :: params
        select case (canonical_model(model_name))
        case ("SYMM_GARCH", "GARCH", "QGARCH", "GJR_SIGNED")
            garch_dist_persist = params%alpha + params%beta
        case ("NAGARCH")
            garch_dist_persist = params%beta + params%alpha*(1.0_dp + params%theta**2)
        case ("GJR_GARCH", "GJR")
            garch_dist_persist = params%alpha + 0.5_dp*params%gamma + params%beta
        case ("EGARCH")
            garch_dist_persist = params%beta
        case default
            garch_dist_persist = params%alpha + params%beta
        end select
    end function garch_dist_persist

    real(dp) function innovation_pdf(dist_id, z, shape, shape2)
        integer, intent(in) :: dist_id
        real(dp), intent(in) :: z, shape, shape2
        select case (dist_id)
        case (dist_normal)
            innovation_pdf = pdf_normal(z)
        case (dist_t)
            innovation_pdf = pdf_t(z, shape)
        case (dist_sech)
            innovation_pdf = pdf_sech(z)
        case (dist_ged)
            innovation_pdf = pdf_ged(z, shape)
        case (dist_laplace)
            innovation_pdf = pdf_laplace(z)
        case (dist_logistic)
            innovation_pdf = pdf_logistic(z)
        case (dist_nig)
            innovation_pdf = pdf_nig(z, shape)
        case (dist_fs_skewt)
            innovation_pdf = pdf_fs_skewt(z, shape, shape2)
        case default
            innovation_pdf = min_pdf
        end select
    end function innovation_pdf

    integer function dist_id_from_name(dist_name)
        character(len=*), intent(in) :: dist_name
        select case (canonical_dist(dist_name))
        case ("NORMAL")
            dist_id_from_name = dist_normal
        case ("T")
            dist_id_from_name = dist_t
        case ("SECH")
            dist_id_from_name = dist_sech
        case ("GED")
            dist_id_from_name = dist_ged
        case ("LAPLACE")
            dist_id_from_name = dist_laplace
        case ("LOGISTIC")
            dist_id_from_name = dist_logistic
        case ("NIG")
            dist_id_from_name = dist_nig
        case ("FS_SKEWT", "FS_SKEW_T", "SKEWT", "SKEW_T")
            dist_id_from_name = dist_fs_skewt
        case default
            error stop "dist_id_from_name: unsupported distribution"
        end select
    end function dist_id_from_name

    character(len=16) function canonical_model(model_name)
        character(len=*), intent(in) :: model_name
        canonical_model = adjustl(model_name)
    end function canonical_model

    character(len=10) function canonical_dist(dist_name)
        character(len=*), intent(in) :: dist_name
        canonical_dist = adjustl(dist_name)
    end function canonical_dist

    real(dp) function sigmoid(x)
        real(dp), intent(in) :: x
        sigmoid = 1.0_dp / (1.0_dp + exp(-clamp(x, -40.0_dp, 40.0_dp)))
    end function sigmoid

    real(dp) function logit_clip(x)
        real(dp), intent(in) :: x
        real(dp) :: xx
        xx = min(max(x, 1.0e-8_dp), 1.0_dp - 1.0e-8_dp)
        logit_clip = log(xx / (1.0_dp - xx))
    end function logit_clip

    real(dp) function atanh_clip(x)
        real(dp), intent(in) :: x
        real(dp) :: xx
        xx = min(max(x, -1.0_dp + 1.0e-8_dp), 1.0_dp - 1.0e-8_dp)
        atanh_clip = 0.5_dp*log((1.0_dp + xx) / (1.0_dp - xx))
    end function atanh_clip

    real(dp) function clamp(x, lo, hi)
        real(dp), intent(in) :: x, lo, hi
        clamp = min(max(x, lo), hi)
    end function clamp

end module garch_fit_dist_mod

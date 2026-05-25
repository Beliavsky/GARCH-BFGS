module garch_sim_mod
    use kind_mod, only: dp
    use garch_types_mod, only: garch_params_t, garch_fit_result_t
    use garch_mod, only: garch_simulate
    use nagarch_mod, only: nagarch_simulate
    use gjr_mod, only: gjr_simulate
    use egarch_mod, only: egarch_simulate
    implicit none
    private

    public :: simulate_garch_model, simulate_garch_fit
    public :: simulate_symm_garch, simulate_nagarch
    public :: simulate_gjr, simulate_gjr_signed, simulate_egarch
    public :: simulate_qgarch, simulate_csgarch, simulate_tgarch
    public :: simulate_symm_garch_pq, simulate_nagarch_pq
    public :: simulate_figarch, simulate_fi_nagarch, simulate_harch
    public :: simulate_avgarch, simulate_aparch
    public :: simulate_midas_hyperbolic, simulate_midas_hyperbolic_asym
    public :: simulate_fgarch_twist

contains

    subroutine simulate_symm_garch(params, n, seed_val, y)
        type(garch_params_t), intent(in) :: params
        integer, intent(in) :: n, seed_val
        real(dp), intent(out) :: y(n)

        call garch_simulate(params%omega, params%alpha, params%beta, &
                            n, seed_val, y)
    end subroutine simulate_symm_garch

    subroutine simulate_nagarch(params, n, seed_val, y)
        type(garch_params_t), intent(in) :: params
        integer, intent(in) :: n, seed_val
        real(dp), intent(out) :: y(n)

        call nagarch_simulate(params%omega, params%alpha, params%beta, &
                              params%theta, n, seed_val, y)
    end subroutine simulate_nagarch

    subroutine simulate_gjr(params, n, seed_val, y)
        type(garch_params_t), intent(in) :: params
        integer, intent(in) :: n, seed_val
        real(dp), intent(out) :: y(n)

        call gjr_simulate(params%omega, params%alpha, params%gamma, &
                          params%beta, n, seed_val, y)
    end subroutine simulate_gjr

    subroutine simulate_gjr_signed(params, n, seed_val, y)
        type(garch_params_t), intent(in) :: params
        integer, intent(in) :: n, seed_val
        real(dp), intent(out) :: y(n)

        call gjr_simulate(params%omega, params%alpha, params%gamma, &
                          params%beta, n, seed_val, y)
    end subroutine simulate_gjr_signed

    subroutine simulate_egarch(params, n, seed_val, y)
        type(garch_params_t), intent(in) :: params
        integer, intent(in) :: n, seed_val
        real(dp), intent(out) :: y(n)

        call egarch_simulate(params%omega, params%alpha, params%gamma, &
                             params%beta, n, seed_val, y)
    end subroutine simulate_egarch

    subroutine simulate_qgarch(params, n, seed_val, y)
        type(garch_params_t), intent(in) :: params
        integer, intent(in) :: n, seed_val
        real(dp), intent(out) :: y(n)
        real(dp) :: h, eps
        integer :: i

        call seed_rng(seed_val)
        h = max((params%omega + params%alpha*params%theta**2) / &
                max(1.0_dp - params%alpha - params%beta, 1.0e-8_dp), 1.0e-12_dp)
        do i = 1, n
            eps = random_normal()
            y(i) = sqrt(h) * eps
            h = max(params%omega + params%alpha*(y(i) - params%theta)**2 + &
                    params%beta*h, 1.0e-12_dp)
        end do
    end subroutine simulate_qgarch

    subroutine simulate_csgarch(params, n, seed_val, y)
        type(garch_params_t), intent(in) :: params
        integer, intent(in) :: n, seed_val
        real(dp), intent(out) :: y(n)
        real(dp) :: backcast, q_prev, h_prev, h, qq, shock
        integer :: i

        call seed_rng(seed_val)
        backcast = max(params%omega / max(1.0_dp - params%extra1, 1.0e-12_dp), 1.0e-12_dp)
        q_prev = backcast
        h_prev = backcast
        shock = backcast

        do i = 1, n
            qq = params%omega + params%extra1*q_prev + params%extra2*(shock - h_prev)
            qq = max(qq, 1.0e-12_dp)
            h = qq + params%alpha*(shock - q_prev) + params%beta*(h_prev - q_prev)
            h = max(h, 1.0e-12_dp)
            y(i) = sqrt(h) * random_normal()
            q_prev = qq
            h_prev = h
            shock = y(i)**2
        end do
    end subroutine simulate_csgarch

    subroutine simulate_tgarch(params, n, seed_val, y)
        type(garch_params_t), intent(in) :: params
        integer, intent(in) :: n, seed_val
        real(dp), intent(out) :: y(n)
        real(dp) :: sigma, eps, news
        integer :: i

        call seed_rng(seed_val)
        sigma = max(params%omega / max(1.0_dp - tgarch_persist(params), 1.0e-8_dp), 1.0e-6_dp)
        do i = 1, n
            eps = random_normal()
            y(i) = sigma * eps
            news = sqrt(max(1.0e-6_dp + eps**2, 1.0e-12_dp)) - params%gamma*eps
            sigma = max(params%omega + params%alpha*news*sigma + params%beta*sigma, 1.0e-6_dp)
        end do
    end subroutine simulate_tgarch

    subroutine simulate_symm_garch_pq(params, n, seed_val, y)
        type(garch_params_t), intent(in) :: params
        integer, intent(in) :: n, seed_val
        real(dp), intent(out) :: y(n)
        real(dp), allocatable :: h(:)
        real(dp) :: backcast, ht, lag_h
        integer :: t, i, lag

        call seed_rng(seed_val)
        allocate(h(n))
        backcast = max(params%omega / max(1.0_dp - sum(params%alpha_lags) - &
                                          sum(params%beta_lags), 1.0e-8_dp), 1.0e-12_dp)
        do t = 1, n
            ht = params%omega
            do i = 1, size(params%alpha_lags)
                lag = t - i
                if (lag >= 1) then
                    ht = ht + params%alpha_lags(i)*y(lag)**2
                else
                    ht = ht + params%alpha_lags(i)*backcast
                end if
            end do
            do i = 1, size(params%beta_lags)
                lag = t - i
                if (lag >= 1) then
                    lag_h = h(lag)
                else
                    lag_h = backcast
                end if
                ht = ht + params%beta_lags(i)*lag_h
            end do
            h(t) = max(ht, 1.0e-12_dp)
            y(t) = sqrt(h(t))*random_normal()
        end do
        deallocate(h)
    end subroutine simulate_symm_garch_pq

    subroutine simulate_nagarch_pq(params, n, seed_val, y)
        type(garch_params_t), intent(in) :: params
        integer, intent(in) :: n, seed_val
        real(dp), intent(out) :: y(n)
        real(dp), allocatable :: h(:)
        real(dp) :: backcast, ht, lag_h, lag_y, sqrth
        integer :: t, i, lag

        call seed_rng(seed_val)
        allocate(h(n))
        backcast = max(params%omega / max(1.0_dp - sum(params%alpha_lags)*(1.0_dp + params%theta**2) - &
                                          sum(params%beta_lags), 1.0e-8_dp), 1.0e-12_dp)
        do t = 1, n
            ht = params%omega
            do i = 1, size(params%alpha_lags)
                lag = t - i
                if (lag >= 1) then
                    lag_h = h(lag)
                    lag_y = y(lag)
                else
                    lag_h = backcast
                    lag_y = 0.0_dp
                end if
                sqrth = sqrt(max(lag_h, 1.0e-12_dp))
                ht = ht + params%alpha_lags(i)*(lag_y - params%theta*sqrth)**2
            end do
            do i = 1, size(params%beta_lags)
                lag = t - i
                if (lag >= 1) then
                    lag_h = h(lag)
                else
                    lag_h = backcast
                end if
                ht = ht + params%beta_lags(i)*lag_h
            end do
            h(t) = max(ht, 1.0e-12_dp)
            y(t) = sqrt(h(t))*random_normal()
        end do
        deallocate(h)
    end subroutine simulate_nagarch_pq

    subroutine simulate_figarch(params, n, seed_val, y)
        type(garch_params_t), intent(in) :: params
        integer, intent(in) :: n, seed_val
        real(dp), intent(out) :: y(n)
        real(dp), allocatable :: lambda(:)
        real(dp) :: backcast, h, omega_tilde, bc_weight
        integer :: t, i, m

        call seed_rng(seed_val)
        m = 1000
        allocate(lambda(m))
        call figarch_weights(params%alpha, params%theta, params%beta, lambda)
        backcast = max(params%omega / max(1.0_dp - sum(lambda), 1.0e-8_dp), 1.0e-12_dp)
        omega_tilde = params%omega / max(1.0_dp - params%beta, 1.0e-8_dp)
        do t = 1, n
            bc_weight = 0.0_dp
            do i = t, m
                bc_weight = bc_weight + lambda(i)
            end do
            h = omega_tilde + bc_weight*backcast
            do i = 1, min(t - 1, m)
                h = h + lambda(i)*y(t - i)**2
            end do
            h = max(h, 1.0e-12_dp)
            y(t) = sqrt(h)*random_normal()
        end do
        deallocate(lambda)
    end subroutine simulate_figarch

    subroutine simulate_fi_nagarch(params, n, seed_val, y)
        type(garch_params_t), intent(in) :: params
        integer, intent(in) :: n, seed_val
        real(dp), intent(out) :: y(n)
        real(dp), allocatable :: lambda(:), news(:)
        real(dp) :: backcast, h, omega_tilde, bc_weight, sqrth, scale
        integer :: t, i, m

        call seed_rng(seed_val)
        m = 1000
        allocate(lambda(m), news(n))
        call figarch_weights(params%alpha, params%theta, params%beta, lambda)
        backcast = max(params%omega / max(1.0_dp - sum(lambda), 1.0e-8_dp), 1.0e-12_dp)
        omega_tilde = params%omega / max(1.0_dp - params%beta, 1.0e-8_dp)
        scale = 1.0_dp + params%twist**2
        do t = 1, n
            bc_weight = 0.0_dp
            do i = t, m
                bc_weight = bc_weight + lambda(i)
            end do
            h = omega_tilde + bc_weight*backcast
            do i = 1, min(t - 1, m)
                h = h + lambda(i)*news(t - i)
            end do
            h = max(h, 1.0e-12_dp)
            sqrth = sqrt(h)
            y(t) = sqrth*random_normal()
            news(t) = (y(t) - params%twist*sqrth)**2 / scale
        end do
        deallocate(lambda, news)
    end subroutine simulate_fi_nagarch

    subroutine simulate_harch(params, n, seed_val, y)
        type(garch_params_t), intent(in) :: params
        integer, intent(in) :: n, seed_val
        real(dp), intent(out) :: y(n)
        real(dp) :: backcast, h
        integer :: t

        call seed_rng(seed_val)
        backcast = max(params%omega / max(1.0_dp - params%alpha - params%gamma - params%beta, 1.0e-8_dp), &
                       1.0e-12_dp)
        do t = 1, n
            h = params%omega + params%alpha*harch_block_y(y, t, 1, backcast) + &
                params%gamma*harch_block_y(y, t, 5, backcast) + &
                params%beta*harch_block_y(y, t, 22, backcast)
            h = max(h, 1.0e-12_dp)
            y(t) = sqrt(h)*random_normal()
        end do
    end subroutine simulate_harch

    subroutine simulate_avgarch(params, n, seed_val, y)
        type(garch_params_t), intent(in) :: params
        integer, intent(in) :: n, seed_val
        real(dp), intent(out) :: y(n)
        real(dp) :: sigma, eps, x, news
        integer :: t

        call seed_rng(seed_val)
        sigma = max(params%omega / max(1.0_dp - avgarch_persist(params), 1.0e-8_dp), 1.0e-6_dp)
        do t = 1, n
            eps = random_normal()
            y(t) = sigma*eps
            x = eps - params%theta
            news = sqrt(max(1.0e-6_dp + x**2, 1.0e-12_dp)) - params%gamma*x
            sigma = max(params%omega + params%alpha*news*sigma + params%beta*sigma, 1.0e-6_dp)
        end do
    end subroutine simulate_avgarch

    subroutine simulate_aparch(params, n, seed_val, y)
        type(garch_params_t), intent(in) :: params
        integer, intent(in) :: n, seed_val
        real(dp), intent(out) :: y(n)
        real(dp) :: sdel, h, eps, term
        integer :: t

        call seed_rng(seed_val)
        sdel = max(params%omega / max(1.0_dp - params%alpha - params%beta, 1.0e-8_dp), 1.0e-12_dp)
        do t = 1, n
            h = max(sdel**(2.0_dp / params%theta), 1.0e-12_dp)
            eps = random_normal()
            y(t) = sqrt(h)*eps
            term = max(abs(y(t)) - params%gamma*y(t), 1.0e-12_dp)**params%theta
            sdel = max(params%omega + params%alpha*term + params%beta*sdel, 1.0e-12_dp)
        end do
    end subroutine simulate_aparch

    subroutine simulate_midas_hyperbolic(params, n, seed_val, y)
        type(garch_params_t), intent(in) :: params
        integer, intent(in) :: n, seed_val
        real(dp), intent(out) :: y(n)
        real(dp) :: weights(22), backcast, x, h
        integer :: t

        call seed_rng(seed_val)
        call midas_weights(params%theta, weights)
        backcast = max(params%omega / max(1.0_dp - params%alpha, 1.0e-8_dp), 1.0e-12_dp)
        do t = 1, n
            x = midas_block_y(y, t, weights, backcast)
            h = max(params%omega + params%alpha*x, 1.0e-12_dp)
            y(t) = sqrt(h)*random_normal()
        end do
    end subroutine simulate_midas_hyperbolic

    subroutine simulate_midas_hyperbolic_asym(params, n, seed_val, y)
        type(garch_params_t), intent(in) :: params
        integer, intent(in) :: n, seed_val
        real(dp), intent(out) :: y(n)
        real(dp) :: weights(22), backcast, x, xneg, h
        integer :: t

        call seed_rng(seed_val)
        call midas_weights(params%theta, weights)
        backcast = max(params%omega / max(1.0_dp - params%alpha - 0.5_dp*params%gamma, 1.0e-8_dp), &
                       1.0e-12_dp)
        do t = 1, n
            call midas_asym_blocks_y(y, t, weights, backcast, x, xneg)
            h = max(params%omega + params%alpha*x + params%gamma*xneg, 1.0e-12_dp)
            y(t) = sqrt(h)*random_normal()
        end do
    end subroutine simulate_midas_hyperbolic_asym

    subroutine simulate_fgarch_twist(params, n, seed_val, y)
        type(garch_params_t), intent(in) :: params
        integer, intent(in) :: n, seed_val
        real(dp), intent(out) :: y(n)
        real(dp) :: h, eps, x, q, moment
        integer :: t

        call seed_rng(seed_val)
        moment = fgarch_twist_moment(params%theta, params%twist)
        h = max(params%omega / max(1.0_dp - params%alpha*moment - params%beta, 1.0e-8_dp), 1.0e-12_dp)
        do t = 1, n
            eps = random_normal()
            y(t) = sqrt(h)*eps
            x = eps - params%theta
            q = abs(x) - params%twist*x
            h = max(params%omega + params%alpha*h*q**2 + params%beta*h, 1.0e-12_dp)
        end do
    end subroutine simulate_fgarch_twist

    subroutine simulate_garch_fit(fit, n, seed_val, y)
        type(garch_fit_result_t), intent(in) :: fit
        integer, intent(in) :: n, seed_val
        real(dp), intent(out) :: y(n)

        call simulate_garch_model(fit%model, fit%params, n, seed_val, y)
    end subroutine simulate_garch_fit

    subroutine simulate_garch_model(model_name, params, n, seed_val, y)
        character(len=*), intent(in) :: model_name
        type(garch_params_t), intent(in) :: params
        integer, intent(in) :: n, seed_val
        real(dp), intent(out) :: y(n)

        select case (trim(adjustl(model_name)))
        case ("SYMM_GARCH", "GARCH")
            call simulate_symm_garch(params, n, seed_val, y)
        case ("NAGARCH")
            call simulate_nagarch(params, n, seed_val, y)
        case ("GJR_GARCH", "GJR")
            call simulate_gjr(params, n, seed_val, y)
        case ("GJR_SIGNED")
            call simulate_gjr_signed(params, n, seed_val, y)
        case ("EGARCH")
            call simulate_egarch(params, n, seed_val, y)
        case ("QGARCH")
            call simulate_qgarch(params, n, seed_val, y)
        case ("CSGARCH")
            call simulate_csgarch(params, n, seed_val, y)
        case ("TGARCH")
            call simulate_tgarch(params, n, seed_val, y)
        case ("SYMM_GARCH_PQ", "SYMM_GARCH21")
            call simulate_symm_garch_pq(params, n, seed_val, y)
        case ("NAGARCH_PQ", "NAGARCH21")
            call simulate_nagarch_pq(params, n, seed_val, y)
        case ("FIGARCH")
            call simulate_figarch(params, n, seed_val, y)
        case ("FI_NAGARCH")
            call simulate_fi_nagarch(params, n, seed_val, y)
        case ("HARCH")
            call simulate_harch(params, n, seed_val, y)
        case ("AVGARCH")
            call simulate_avgarch(params, n, seed_val, y)
        case ("APARCH")
            call simulate_aparch(params, n, seed_val, y)
        case ("MIDAS_HYPER", "MIDAS_HYPERBOLIC")
            call simulate_midas_hyperbolic(params, n, seed_val, y)
        case ("MIDAS_ASYM", "MIDAS_HYPERBOLIC_ASYM")
            call simulate_midas_hyperbolic_asym(params, n, seed_val, y)
        case ("FGARCH_TWIST")
            call simulate_fgarch_twist(params, n, seed_val, y)
        case default
            error stop "garch_sim_mod: no simulator for requested model"
        end select
    end subroutine simulate_garch_model

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

    real(dp) function tgarch_persist(params)
        type(garch_params_t), intent(in) :: params
        real(dp) :: kappa

        kappa = sqrt(2.0_dp / acos(-1.0_dp))
        tgarch_persist = params%alpha*kappa + params%beta
    end function tgarch_persist

    real(dp) function avgarch_persist(params)
        type(garch_params_t), intent(in) :: params

        avgarch_persist = params%alpha*avgarch_kappa(params%gamma, params%theta) + params%beta
    end function avgarch_persist

    real(dp) function avgarch_kappa(eta, loc)
        real(dp), intent(in) :: eta, loc
        real(dp) :: phi, cdf, pi

        pi = acos(-1.0_dp)
        phi = exp(-0.5_dp*loc**2) / sqrt(2.0_dp*pi)
        cdf = 0.5_dp * (1.0_dp + erf(loc / sqrt(2.0_dp)))
        avgarch_kappa = 2.0_dp*phi + loc*(2.0_dp*cdf - 1.0_dp) + eta*loc
    end function avgarch_kappa

    subroutine figarch_weights(phi, d, beta, lambda)
        real(dp), intent(in) :: phi, d, beta
        real(dp), intent(out) :: lambda(:)
        real(dp) :: delta_prev, delta_cur
        integer :: i

        if (size(lambda) < 1) return
        delta_prev = d
        lambda(1) = d - phi + beta
        do i = 2, size(lambda)
            delta_cur = ((real(i, dp) - 1.0_dp - d) / real(i, dp))*delta_prev
            lambda(i) = beta*lambda(i - 1) + delta_cur - phi*delta_prev
            delta_prev = delta_cur
        end do
    end subroutine figarch_weights

    real(dp) function harch_block_y(y, t, lag, backcast)
        real(dp), intent(in) :: y(:), backcast
        integer, intent(in) :: t, lag
        integer :: j, idx

        harch_block_y = 0.0_dp
        do j = 1, lag
            idx = t - j
            if (idx >= 1) then
                harch_block_y = harch_block_y + y(idx)**2
            else
                harch_block_y = harch_block_y + backcast
            end if
        end do
        harch_block_y = harch_block_y / real(lag, dp)
    end function harch_block_y

    subroutine midas_weights(theta, weights)
        real(dp), intent(in) :: theta
        real(dp), intent(out) :: weights(:)
        real(dp) :: raw(size(weights)), ratio, sum_raw
        integer :: i

        raw(1) = theta
        do i = 2, size(weights)
            ratio = (real(i - 1, dp) + theta) / real(i, dp)
            raw(i) = raw(i - 1)*ratio
        end do
        sum_raw = sum(raw)
        weights = raw / max(sum_raw, 1.0e-12_dp)
    end subroutine midas_weights

    real(dp) function midas_block_y(y, t, weights, backcast)
        real(dp), intent(in) :: y(:), weights(:), backcast
        integer, intent(in) :: t
        integer :: i, idx

        midas_block_y = 0.0_dp
        do i = 1, size(weights)
            idx = t - i
            if (idx >= 1) then
                midas_block_y = midas_block_y + weights(i)*y(idx)**2
            else
                midas_block_y = midas_block_y + weights(i)*backcast
            end if
        end do
    end function midas_block_y

    subroutine midas_asym_blocks_y(y, t, weights, backcast, x, xneg)
        real(dp), intent(in) :: y(:), weights(:), backcast
        integer, intent(in) :: t
        real(dp), intent(out) :: x, xneg
        real(dp) :: lag_sq, ind
        integer :: i, idx

        x = 0.0_dp
        xneg = 0.0_dp
        do i = 1, size(weights)
            idx = t - i
            if (idx >= 1) then
                lag_sq = y(idx)**2
                ind = merge(1.0_dp, 0.0_dp, y(idx) < 0.0_dp)
            else
                lag_sq = backcast
                ind = 0.5_dp
            end if
            x = x + weights(i)*lag_sq
            xneg = xneg + weights(i)*ind*lag_sq
        end do
    end subroutine midas_asym_blocks_y

    real(dp) function fgarch_twist_moment(theta, c)
        real(dp), intent(in) :: theta, c
        real(dp) :: phi, ph, a, b, one_minus_c, one_plus_c

        phi = 0.5_dp * (1.0_dp + erf(theta / sqrt(2.0_dp)))
        ph = exp(-0.5_dp*theta**2) / sqrt(2.0_dp*acos(-1.0_dp))
        a = (1.0_dp + theta**2)*(1.0_dp - phi) - theta*ph
        b = (1.0_dp + theta**2)*phi + theta*ph
        one_minus_c = 1.0_dp - c
        one_plus_c = 1.0_dp + c
        fgarch_twist_moment = one_minus_c**2*a + one_plus_c**2*b
    end function fgarch_twist_moment

end module garch_sim_mod

module garch_fit_mod
    use kind_mod,       only: dp
    use math_const_mod, only: log_sqrt_2pi, sqrt2
    use garch_types_mod, only: garch_params_t
    use garch_mod,      only: garch_set_data, garch_obj, garch_transform, garch_inv_transform
    use nagarch_mod,    only: nagarch_set_data, nagarch_obj, nagarch_transform, nagarch_inv_transform
    use gjr_mod,        only: gjr_set_data, gjr_obj, gjr_transform, gjr_inv_transform, &
                              gjr_signed_obj, gjr_signed_transform, gjr_signed_inv_transform
    use egarch_mod,     only: egarch_set_data, egarch_obj, egarch_transform, egarch_inv_transform
    use fgarch_mod,     only: fg_dist_normal, fgarch_set_data, fgarch_set_dist, &
                              fgarch_twist_obj, fgarch_twist_transform, fgarch_twist_inv_transform, &
                              fgarch_twist_vol_ann, fgarch_twist_skew_kurt
    use bfgs_mod,       only: bfgs_minimize
    implicit none
    private

    public :: fit_symm_garch, fit_qgarch, fit_figarch, fit_fi_nagarch, fit_nagarch, fit_rgarch, fit_carr_park
    public :: fit_regarch1, fit_regarch2, fit_rgarch_meas
    public :: fit_nagarch_range, fit_gjr, fit_gjr_signed, fit_egarch, fit_aparch, fit_harch
    public :: fit_riskmetrics2006, fit_midas_hyperbolic, fit_midas_hyperbolic_asym
    public :: fit_fgarch_twist, fit_fgarch_twist_range, fit_ewma
    public :: fit_aewma_nag, fit_aewma_twist
    public :: garch_skew_kurt, qgarch_skew_kurt, nagarch_skew_kurt, rgarch_skew_kurt, carr_park_skew_kurt
    public :: regarch1_skew_kurt, regarch2_skew_kurt
    public :: rgarch_meas_skew_kurt
    public :: nagarch_range_skew_kurt
    public :: gjr_skew_kurt, egarch_skew_kurt, aparch_skew_kurt, harch_skew_kurt
    public :: figarch_skew_kurt, fi_nagarch_skew_kurt
    public :: riskmetrics2006_skew_kurt, midas_hyperbolic_skew_kurt, midas_hyperbolic_asym_skew_kurt
    public :: ewma_skew_kurt
    public :: aewma_nag_skew_kurt, aewma_twist_skew_kurt, fgarch_twist_range_skew_kurt
    public :: symm_garch_persist, qgarch_persist, figarch_persist, fi_nagarch_persist, nagarch_persist
    public :: gjr_persist, egarch_persist, aparch_persist, harch_persist
    public :: riskmetrics2006_persist, midas_hyperbolic_persist, midas_hyperbolic_asym_persist
    public :: aparch_mean_variance, qgarch_mean_variance, figarch_variance, fi_nagarch_variance
    public :: riskmetrics2006_variance
    public :: fgarch_twist_moment, fgarch_twist_persist, ewma_persist
    public :: aewma_nag_persist, aewma_twist_persist, rgarch_persist, carr_park_persist
    public :: regarch1_persist, regarch2_persist
    public :: rgarch_meas_persist

    integer, parameter :: n_start = 4
    integer, parameter :: ewma_np = 1
    integer, parameter :: aewma_nag_np = 3
    integer, parameter :: aewma_twist_np = 4
    integer, parameter :: symm_garch_np = 3
    integer, parameter :: qgarch_np = 4
    integer, parameter :: figarch_np = 4
    integer, parameter :: fi_nagarch_np = 5
    integer, parameter :: figarch_trunc_lag = 1000
    integer, parameter :: nagarch_np = 4
    integer, parameter :: rgarch_np = 3
    integer, parameter :: carr_park_np = 3
    integer, parameter :: regarch1_np = 4
    integer, parameter :: regarch2_np = 7
    integer, parameter :: rgarch_meas_np = 9
    integer, parameter :: nagarch_range_np = 5
    integer, parameter :: gjr_np = 4
    integer, parameter :: egarch_np = 4
    integer, parameter :: aparch_np = 5
    integer, parameter :: harch_np = 4
    integer, parameter :: midas_hyperbolic_np = 3
    integer, parameter :: midas_hyperbolic_asym_np = 4
    integer, parameter :: midas_hyperbolic_m = 22
    integer, parameter :: rm2006_kmax = 14
    integer, parameter :: fgarch_twist_np = 5
    integer, parameter :: fgarch_twist_range_np = 6
    real(dp), allocatable, save :: ewma_obs(:)
    real(dp), allocatable, save :: qgarch_obs(:)
    real(dp), allocatable, save :: figarch_obs(:)
    real(dp), allocatable, save :: aparch_obs(:)
    real(dp), allocatable, save :: harch_obs(:)
    real(dp), allocatable, save :: midas_hyperbolic_obs(:)
    real(dp), allocatable, save :: nag_range_obs(:), nag_range_x(:)
    real(dp), allocatable, save :: rgarch_obs(:), rgarch_x(:)
    real(dp), allocatable, save :: carr_park_x(:)
    real(dp), allocatable, save :: regarch_obs(:), regarch_log_range(:)
    integer, save :: ewma_nobs = 0
    integer, save :: qgarch_nobs = 0
    integer, save :: figarch_nobs = 0
    integer, save :: aparch_nobs = 0
    integer, save :: harch_nobs = 0
    integer, save :: midas_hyperbolic_nobs = 0
    integer, save :: nag_range_nobs = 0
    integer, save :: rgarch_nobs = 0
    integer, save :: carr_park_nobs = 0
    integer, save :: regarch_nobs = 0
    real(dp), parameter :: regarch_log_range_mean = 0.43_dp
    real(dp), parameter :: regarch_log_range_sd = 0.29_dp

contains

    subroutine fit_ewma(y, max_iter, gtol, f_best, params, niter_best, converged_best)
        real(dp), intent(in)  :: y(:), gtol
        integer,  intent(in)  :: max_iter
        real(dp), intent(out) :: f_best
        type(garch_params_t), intent(out) :: params
        integer,  intent(out) :: niter_best
        logical,  intent(out) :: converged_best
        real(dp), parameter :: starts(n_start) = [0.90_dp, 0.94_dp, 0.97_dp, 0.99_dp]
        real(dp) :: p(ewma_np), p0(ewma_np), p_best(ewma_np), f_try
        integer :: istart, niter_try
        logical :: converged_try

        call ewma_set_data(y)
        f_best = huge(1.0_dp)
        p_best = 0.0_dp
        niter_best = 0
        converged_best = .false.

        do istart = 1, n_start
            call ewma_inv_transform(starts(istart), p0)
            p = p0
            call bfgs_minimize(ewma_obj, p, ewma_np, max_iter, gtol, f_try, niter_try, converged_try)
            if (f_try < f_best) then
                f_best = f_try
                p_best = p
                niter_best = niter_try
                converged_best = converged_try
            end if
        end do

        params = garch_params_t()
        call ewma_transform(p_best, params%beta)
        params%omega = 0.0_dp
        params%alpha = 1.0_dp - params%beta
    end subroutine fit_ewma

    subroutine fit_aewma_nag(y, max_iter, gtol, f_best, params, niter_best, converged_best)
        real(dp), intent(in)  :: y(:), gtol
        integer,  intent(in)  :: max_iter
        real(dp), intent(out) :: f_best
        type(garch_params_t), intent(out) :: params
        integer,  intent(out) :: niter_best
        logical,  intent(out) :: converged_best
        real(dp), parameter :: start_lambda(n_start) = [0.90_dp, 0.94_dp, 0.97_dp, 0.99_dp]
        real(dp), parameter :: start_theta(n_start) = [0.0_dp, 0.5_dp, 1.0_dp, -0.5_dp]
        real(dp) :: p(aewma_nag_np), p0(aewma_nag_np), p_best(aewma_nag_np), f_try
        real(dp) :: scale0
        integer :: istart, niter_try
        logical :: converged_try

        call ewma_set_data(y)
        f_best = huge(1.0_dp)
        p_best = 0.0_dp
        niter_best = 0
        converged_best = .false.

        do istart = 1, n_start
            scale0 = 1.0_dp / (1.0_dp + start_theta(istart)**2)
            call aewma_nag_inv_transform(start_lambda(istart), start_theta(istart), scale0, p0)
            p = p0
            call bfgs_minimize(aewma_nag_obj, p, aewma_nag_np, max_iter, gtol, f_try, niter_try, converged_try)
            if (f_try < f_best) then
                f_best = f_try
                p_best = p
                niter_best = niter_try
                converged_best = converged_try
            end if
        end do

        params = garch_params_t()
        call aewma_nag_transform(p_best, params%beta, params%theta, params%scale)
        params%omega = 0.0_dp
        params%alpha = 1.0_dp - params%beta
    end subroutine fit_aewma_nag

    subroutine fit_aewma_twist(y, max_iter, gtol, f_best, params, niter_best, converged_best)
        real(dp), intent(in)  :: y(:), gtol
        integer,  intent(in)  :: max_iter
        real(dp), intent(out) :: f_best
        type(garch_params_t), intent(out) :: params
        integer,  intent(out) :: niter_best
        logical,  intent(out) :: converged_best
        real(dp), parameter :: start_lambda(n_start) = [0.90_dp, 0.94_dp, 0.97_dp, 0.99_dp]
        real(dp), parameter :: start_theta(n_start) = [0.0_dp, 0.5_dp, 1.0_dp, 1.0_dp]
        real(dp), parameter :: start_twist(n_start) = [0.0_dp, 0.0_dp, -0.25_dp, -0.50_dp]
        real(dp) :: p(aewma_twist_np), p0(aewma_twist_np), p_best(aewma_twist_np), f_try
        real(dp) :: scale0
        integer :: istart, niter_try
        logical :: converged_try

        call ewma_set_data(y)
        f_best = huge(1.0_dp)
        p_best = 0.0_dp
        niter_best = 0
        converged_best = .false.

        do istart = 1, n_start
            scale0 = 1.0_dp / fgarch_twist_moment(start_theta(istart), start_twist(istart))
            call aewma_twist_inv_transform(start_lambda(istart), start_theta(istart), &
                                           start_twist(istart), scale0, p0)
            p = p0
            call bfgs_minimize(aewma_twist_obj, p, aewma_twist_np, max_iter, gtol, f_try, niter_try, converged_try)
            if (f_try < f_best) then
                f_best = f_try
                p_best = p
                niter_best = niter_try
                converged_best = converged_try
            end if
        end do

        params = garch_params_t()
        call aewma_twist_transform(p_best, params%beta, params%theta, params%twist, params%scale)
        params%omega = 0.0_dp
        params%alpha = 1.0_dp - params%beta
    end subroutine fit_aewma_twist

    subroutine fit_symm_garch(y, max_iter, gtol, f_best, params, niter_best, converged_best)
        real(dp), intent(in)  :: y(:), gtol
        integer,  intent(in)  :: max_iter
        real(dp), intent(out) :: f_best
        type(garch_params_t), intent(out) :: params
        integer,  intent(out) :: niter_best
        logical,  intent(out) :: converged_best
        real(dp), parameter :: starts(symm_garch_np,n_start) = reshape( &
            [1.0e-6_dp, 0.04_dp, 0.90_dp, &
             1.0e-6_dp, 0.06_dp, 0.88_dp, &
             5.0e-6_dp, 0.08_dp, 0.85_dp, &
             5.0e-6_dp, 0.10_dp, 0.80_dp], [symm_garch_np,n_start])
        real(dp) :: p(symm_garch_np), p0(symm_garch_np), p_best(symm_garch_np), f_try
        integer :: istart, niter_try, n
        logical :: converged_try

        n = size(y)
        call garch_set_data(y, n)
        f_best = huge(1.0_dp)
        p_best = 0.0_dp
        niter_best = 0
        converged_best = .false.

        do istart = 1, n_start
            call garch_inv_transform(starts(1,istart), starts(2,istart), starts(3,istart), p0)
            p = p0
            call bfgs_minimize(garch_obj, p, symm_garch_np, max_iter, gtol, f_try, niter_try, converged_try)
            if (f_try < f_best) then
                f_best = f_try
                p_best = p
                niter_best = niter_try
                converged_best = converged_try
            end if
        end do

        params = garch_params_t()
        call garch_transform(p_best, params%omega, params%alpha, params%beta)
    end subroutine fit_symm_garch

    subroutine fit_qgarch(y, max_iter, gtol, f_best, params, niter_best, converged_best)
        real(dp), intent(in)  :: y(:), gtol
        integer,  intent(in)  :: max_iter
        real(dp), intent(out) :: f_best
        type(garch_params_t), intent(out) :: params
        integer,  intent(out) :: niter_best
        logical,  intent(out) :: converged_best
        real(dp), parameter :: start_alpha(n_start) = [0.05_dp, 0.08_dp, 0.12_dp, 0.04_dp]
        real(dp), parameter :: start_beta(n_start)  = [0.90_dp, 0.86_dp, 0.80_dp, 0.94_dp]
        real(dp), parameter :: start_shift(n_start) = [0.0_dp, -0.005_dp, 0.005_dp, -0.010_dp]
        real(dp) :: p(qgarch_np), p0(qgarch_np), p_best(qgarch_np), f_try
        real(dp) :: var_y, omega0
        integer :: istart, niter_try
        logical :: converged_try

        call qgarch_set_data(y)
        var_y = max(sum(y**2) / real(size(y), dp), 1.0e-12_dp)
        f_best = huge(1.0_dp)
        p_best = 0.0_dp
        niter_best = 0
        converged_best = .false.

        if (qgarch_nagarch_start(y, max_iter, gtol, p0)) then
            p = p0
            call bfgs_minimize(qgarch_obj, p, qgarch_np, max_iter, gtol, f_try, niter_try, converged_try)
            if (f_try < f_best) then
                f_best = f_try
                p_best = p
                niter_best = niter_try
                converged_best = converged_try
            end if
            if (converged_try) then
                params = garch_params_t()
                call qgarch_transform(p_best, params%omega, params%alpha, params%beta, params%theta)
                return
            end if
        end if

        do istart = 1, n_start
            omega0 = max((1.0_dp - start_alpha(istart) - start_beta(istart))*var_y - &
                         start_alpha(istart)*start_shift(istart)**2, 1.0e-12_dp)
            call qgarch_inv_transform(omega0, start_alpha(istart), start_beta(istart), start_shift(istart), p0)
            p = p0
            call bfgs_minimize(qgarch_obj, p, qgarch_np, max_iter, gtol, f_try, niter_try, converged_try)
            if (f_try < f_best) then
                f_best = f_try
                p_best = p
                niter_best = niter_try
                converged_best = converged_try
            end if
        end do

        params = garch_params_t()
        call qgarch_transform(p_best, params%omega, params%alpha, params%beta, params%theta)
    end subroutine fit_qgarch

    subroutine fit_figarch(y, max_iter, gtol, f_best, params, niter_best, converged_best)
        real(dp), intent(in)  :: y(:), gtol
        integer,  intent(in)  :: max_iter
        real(dp), intent(out) :: f_best
        type(garch_params_t), intent(out) :: params
        integer,  intent(out) :: niter_best
        logical,  intent(out) :: converged_best
        real(dp), parameter :: start_phi(n_start)  = [0.10_dp, 0.20_dp, 0.05_dp, 0.30_dp]
        real(dp), parameter :: start_d(n_start)    = [0.40_dp, 0.60_dp, 0.20_dp, 0.80_dp]
        real(dp), parameter :: start_beta(n_start) = [0.40_dp, 0.50_dp, 0.20_dp, 0.60_dp]
        real(dp) :: p(figarch_np), p0(figarch_np), p_best(figarch_np), f_try
        real(dp) :: var_y, omega0, beta0
        integer :: istart, niter_try
        logical :: converged_try

        call figarch_set_data(y)
        var_y = max(sum(y**2) / real(size(y), dp), 1.0e-12_dp)
        f_best = huge(1.0_dp)
        p_best = 0.0_dp
        niter_best = 0
        converged_best = .false.

        do istart = 1, n_start
            beta0 = min(start_beta(istart), 0.95_dp*(start_d(istart) + start_phi(istart)))
            omega0 = max((1.0_dp - beta0)*0.10_dp*var_y, 1.0e-12_dp)
            call figarch_inv_transform(omega0, start_phi(istart), start_d(istart), beta0, p0)
            p = p0
            call bfgs_minimize(figarch_obj, p, figarch_np, max_iter, gtol, f_try, niter_try, converged_try)
            if (f_try < f_best) then
                f_best = f_try
                p_best = p
                niter_best = niter_try
                converged_best = converged_try
            end if
        end do

        params = garch_params_t()
        call figarch_transform(p_best, params%omega, params%alpha, params%theta, params%beta)
    end subroutine fit_figarch

    subroutine fit_fi_nagarch(y, max_iter, gtol, f_best, params, niter_best, converged_best)
        real(dp), intent(in)  :: y(:), gtol
        integer,  intent(in)  :: max_iter
        real(dp), intent(out) :: f_best
        type(garch_params_t), intent(out) :: params
        integer,  intent(out) :: niter_best
        logical,  intent(out) :: converged_best
        real(dp), parameter :: start_shift(n_start) = [0.0_dp, 0.5_dp, -0.5_dp, 1.0_dp]
        type(garch_params_t) :: fig_params
        real(dp) :: p(fi_nagarch_np), p0(fi_nagarch_np), p_best(fi_nagarch_np), f_try, f_fig
        integer :: istart, niter_try, niter_fig
        logical :: converged_try, converged_fig

        call fit_figarch(y, max_iter, gtol, f_fig, fig_params, niter_fig, converged_fig)
        call figarch_set_data(y)
        f_best = huge(1.0_dp)
        p_best = 0.0_dp
        niter_best = 0
        converged_best = .false.

        do istart = 1, n_start
            call fi_nagarch_inv_transform(fig_params%omega, fig_params%alpha, fig_params%theta, &
                                          fig_params%beta, start_shift(istart), p0)
            p = p0
            call bfgs_minimize(fi_nagarch_obj, p, fi_nagarch_np, max_iter, gtol, f_try, niter_try, converged_try)
            if (f_try < f_best) then
                f_best = f_try
                p_best = p
                niter_best = niter_try + niter_fig
                converged_best = converged_try
            end if
        end do

        params = garch_params_t()
        call fi_nagarch_transform(p_best, params%omega, params%alpha, params%theta, params%beta, params%twist)
    end subroutine fit_fi_nagarch

    subroutine fit_harch(y, max_iter, gtol, f_best, params, niter_best, converged_best)
        real(dp), intent(in)  :: y(:), gtol
        integer,  intent(in)  :: max_iter
        real(dp), intent(out) :: f_best
        type(garch_params_t), intent(out) :: params
        integer,  intent(out) :: niter_best
        logical,  intent(out) :: converged_best
        real(dp), parameter :: starts(harch_np,n_start) = reshape( &
            [1.0e-6_dp, 0.10_dp, 0.40_dp, 0.40_dp, &
             1.0e-6_dp, 0.00_dp, 0.45_dp, 0.40_dp, &
             1.0e-6_dp, 0.10_dp, 0.20_dp, 0.60_dp, &
             5.0e-6_dp, 0.05_dp, 0.30_dp, 0.50_dp], [harch_np,n_start])
        real(dp) :: p(harch_np), p0(harch_np), p_best(harch_np), f_try, omega0
        integer :: istart, niter_try, n
        logical :: converged_try

        n = size(y)
        call harch_set_data(y)
        f_best = huge(1.0_dp)
        p_best = 0.0_dp
        niter_best = 0
        converged_best = .false.

        do istart = 1, n_start
            omega0 = max((1.0_dp - starts(2,istart) - starts(3,istart) - starts(4,istart)) * &
                         sum(y**2) / real(n, dp), 1.0e-12_dp)
            call harch_inv_transform(omega0, starts(2,istart), starts(3,istart), starts(4,istart), p0)
            p = p0
            call bfgs_minimize(harch_obj, p, harch_np, max_iter, gtol, f_try, niter_try, converged_try)
            if (f_try < f_best) then
                f_best = f_try
                p_best = p
                niter_best = niter_try
                converged_best = converged_try
            end if
        end do

        params = garch_params_t()
        call harch_transform(p_best, params%omega, params%alpha, params%gamma, params%beta)
    end subroutine fit_harch

    subroutine fit_riskmetrics2006(y, f_best, params, niter_best, converged_best)
        real(dp), intent(in)  :: y(:)
        real(dp), intent(out) :: f_best
        type(garch_params_t), intent(out) :: params
        integer, intent(out) :: niter_best
        logical, intent(out) :: converged_best
        real(dp), allocatable :: variance(:)
        integer :: t, n

        n = size(y)
        allocate(variance(n))
        call riskmetrics2006_variance(y, variance)
        f_best = real(n, dp) * log_sqrt_2pi
        do t = 1, n
            f_best = f_best + 0.5_dp * (log(max(variance(t), 1.0e-12_dp)) + y(t)**2 / max(variance(t), 1.0e-12_dp))
        end do
        f_best = f_best / real(n, dp)
        params = garch_params_t()
        niter_best = 0
        converged_best = .true.
        deallocate(variance)
    end subroutine fit_riskmetrics2006

    subroutine fit_midas_hyperbolic(y, max_iter, gtol, f_best, params, niter_best, converged_best)
        real(dp), intent(in)  :: y(:), gtol
        integer,  intent(in)  :: max_iter
        real(dp), intent(out) :: f_best
        type(garch_params_t), intent(out) :: params
        integer,  intent(out) :: niter_best
        logical,  intent(out) :: converged_best
        real(dp), parameter :: start_alpha(n_start) = [0.80_dp, 0.90_dp, 0.95_dp, 0.98_dp]
        real(dp), parameter :: start_theta(n_start) = [0.10_dp, 0.50_dp, 0.80_dp, 0.90_dp]
        real(dp) :: p(midas_hyperbolic_np), p0(midas_hyperbolic_np), p_best(midas_hyperbolic_np), f_try
        real(dp) :: omega0, var_y
        integer :: istart, niter_try
        logical :: converged_try

        call midas_hyperbolic_set_data(y)
        var_y = max(sum(y**2) / real(size(y), dp), 1.0e-12_dp)
        f_best = huge(1.0_dp)
        p_best = 0.0_dp
        niter_best = 0
        converged_best = .false.

        do istart = 1, n_start
            omega0 = max((1.0_dp - min(start_alpha(istart), 0.99_dp))*var_y, 1.0e-12_dp)
            call midas_hyperbolic_inv_transform(omega0, start_alpha(istart), start_theta(istart), p0)
            p = p0
            call bfgs_minimize(midas_hyperbolic_obj, p, midas_hyperbolic_np, max_iter, gtol, &
                               f_try, niter_try, converged_try)
            if (f_try < f_best) then
                f_best = f_try
                p_best = p
                niter_best = niter_try
                converged_best = converged_try
            end if
        end do

        params = garch_params_t()
        call midas_hyperbolic_transform(p_best, params%omega, params%alpha, params%theta)
    end subroutine fit_midas_hyperbolic

    subroutine fit_midas_hyperbolic_asym(y, max_iter, gtol, f_best, params, niter_best, converged_best)
        real(dp), intent(in)  :: y(:), gtol
        integer,  intent(in)  :: max_iter
        real(dp), intent(out) :: f_best
        type(garch_params_t), intent(out) :: params
        integer,  intent(out) :: niter_best
        logical,  intent(out) :: converged_best
        real(dp), parameter :: start_alpha(4) = [0.80_dp, 0.90_dp, 0.95_dp, 0.98_dp]
        real(dp), parameter :: start_theta(4) = [0.10_dp, 0.50_dp, 0.80_dp, 0.90_dp]
        real(dp), parameter :: start_gamma(3) = [0.0_dp, 0.5_dp, 0.9_dp]
        real(dp) :: p(midas_hyperbolic_asym_np), p0(midas_hyperbolic_asym_np), p_best(midas_hyperbolic_asym_np)
        real(dp) :: f_try, omega0, var_y, total
        integer :: ia, itheta, ig, niter_try
        logical :: converged_try

        call midas_hyperbolic_set_data(y)
        var_y = max(sum(y**2) / real(size(y), dp), 1.0e-12_dp)
        f_best = huge(1.0_dp)
        p_best = 0.0_dp
        niter_best = 0
        converged_best = .false.

        do ia = 1, size(start_alpha)
            do itheta = 1, size(start_theta)
                do ig = 1, size(start_gamma)
                    total = start_alpha(ia) + 0.5_dp*start_gamma(ig)
                    omega0 = max((1.0_dp - min(total, 0.99_dp))*var_y, 1.0e-12_dp)
                    call midas_hyperbolic_asym_inv_transform(omega0, start_alpha(ia), start_gamma(ig), &
                                                             start_theta(itheta), p0)
                    p = p0
                    call bfgs_minimize(midas_hyperbolic_asym_obj, p, midas_hyperbolic_asym_np, max_iter, gtol, &
                                       f_try, niter_try, converged_try)
                    if (f_try < f_best) then
                        f_best = f_try
                        p_best = p
                        niter_best = niter_try
                        converged_best = converged_try
                    end if
                end do
            end do
        end do

        params = garch_params_t()
        call midas_hyperbolic_asym_transform(p_best, params%omega, params%alpha, params%gamma, params%theta)
    end subroutine fit_midas_hyperbolic_asym

    subroutine fit_nagarch(y, max_iter, gtol, f_best, params, niter_best, converged_best)
        real(dp), intent(in)  :: y(:), gtol
        integer,  intent(in)  :: max_iter
        real(dp), intent(out) :: f_best
        type(garch_params_t), intent(out) :: params
        integer,  intent(out) :: niter_best
        logical,  intent(out) :: converged_best
        real(dp), parameter :: starts(nagarch_np,n_start) = reshape( &
            [1.0e-5_dp, 0.04_dp, 0.90_dp, -0.50_dp, &
             1.0e-5_dp, 0.06_dp, 0.88_dp,  0.00_dp, &
             1.0e-5_dp, 0.06_dp, 0.88_dp,  0.50_dp, &
             5.0e-6_dp, 0.10_dp, 0.80_dp,  0.50_dp], [nagarch_np,n_start])
        real(dp) :: p(nagarch_np), p0(nagarch_np), p_best(nagarch_np), f_try
        integer :: istart, niter_try, n
        logical :: converged_try

        n = size(y)
        call nagarch_set_data(y, n)
        f_best = huge(1.0_dp)
        p_best = 0.0_dp
        niter_best = 0
        converged_best = .false.

        do istart = 1, n_start
            call nagarch_inv_transform(starts(1,istart), starts(2,istart), &
                                       starts(3,istart), starts(4,istart), p0)
            p = p0
            call bfgs_minimize(nagarch_obj, p, nagarch_np, max_iter, gtol, f_try, niter_try, converged_try)
            if (f_try < f_best) then
                f_best = f_try
                p_best = p
                niter_best = niter_try
                converged_best = converged_try
            end if
        end do

        params = garch_params_t()
        call nagarch_transform(p_best, params%omega, params%alpha, params%beta, params%theta)
    end subroutine fit_nagarch

    subroutine fit_rgarch(y, range_var, max_iter, gtol, f_best, params, niter_best, converged_best)
        real(dp), intent(in)  :: y(:), range_var(:), gtol
        integer,  intent(in)  :: max_iter
        real(dp), intent(out) :: f_best
        type(garch_params_t), intent(out) :: params
        integer,  intent(out) :: niter_best
        logical,  intent(out) :: converged_best
        real(dp), parameter :: starts(rgarch_np,n_start) = reshape( &
            [1.0e-6_dp, 0.10_dp, 0.85_dp, &
             1.0e-6_dp, 0.20_dp, 0.70_dp, &
             5.0e-6_dp, 0.30_dp, 0.60_dp, &
             5.0e-6_dp, 0.50_dp, 0.40_dp], [rgarch_np,n_start])
        real(dp) :: p(rgarch_np), p0(rgarch_np), p_best(rgarch_np), f_try
        integer :: istart, niter_try
        logical :: converged_try

        call rgarch_set_data(y, range_var)
        f_best = huge(1.0_dp)
        p_best = 0.0_dp
        niter_best = 0
        converged_best = .false.

        do istart = 1, n_start
            call rgarch_inv_transform(starts(1,istart), starts(2,istart), starts(3,istart), p0)
            p = p0
            call bfgs_minimize(rgarch_obj, p, rgarch_np, max_iter, gtol, f_try, niter_try, converged_try)
            if (f_try < f_best) then
                f_best = f_try
                p_best = p
                niter_best = niter_try
                converged_best = converged_try
            end if
        end do

        params = garch_params_t()
        call rgarch_transform(p_best, params%omega, params%alpha, params%beta)
    end subroutine fit_rgarch

    subroutine fit_carr_park(range_var, max_iter, gtol, f_best, params, niter_best, converged_best)
        real(dp), intent(in)  :: range_var(:), gtol
        integer,  intent(in)  :: max_iter
        real(dp), intent(out) :: f_best
        type(garch_params_t), intent(out) :: params
        integer,  intent(out) :: niter_best
        logical,  intent(out) :: converged_best
        real(dp), parameter :: starts(carr_park_np,n_start) = reshape( &
            [1.0e-6_dp, 0.10_dp, 0.85_dp, &
             1.0e-6_dp, 0.20_dp, 0.70_dp, &
             5.0e-6_dp, 0.30_dp, 0.60_dp, &
             5.0e-6_dp, 0.50_dp, 0.40_dp], [carr_park_np,n_start])
        real(dp) :: p(carr_park_np), p0(carr_park_np), p_best(carr_park_np), f_try
        integer :: istart, niter_try
        logical :: converged_try

        call carr_park_set_data(range_var)
        f_best = huge(1.0_dp)
        p_best = 0.0_dp
        niter_best = 0
        converged_best = .false.

        do istart = 1, n_start
            call rgarch_inv_transform(starts(1,istart), starts(2,istart), starts(3,istart), p0)
            p = p0
            call bfgs_minimize(carr_park_obj, p, carr_park_np, max_iter, gtol, &
                               f_try, niter_try, converged_try)
            if (f_try < f_best) then
                f_best = f_try
                p_best = p
                niter_best = niter_try
                converged_best = converged_try
            end if
        end do

        params = garch_params_t()
        call rgarch_transform(p_best, params%omega, params%alpha, params%beta)
    end subroutine fit_carr_park

    subroutine fit_regarch1(y, log_range, max_iter, gtol, f_best, params, niter_best, converged_best)
        real(dp), intent(in)  :: y(:), log_range(:), gtol
        integer,  intent(in)  :: max_iter
        real(dp), intent(out) :: f_best
        type(garch_params_t), intent(out) :: params
        integer,  intent(out) :: niter_best
        logical,  intent(out) :: converged_best
        real(dp), parameter :: start_alpha(n_start) = [0.20_dp, 0.10_dp, 0.30_dp, 0.10_dp]
        real(dp), parameter :: start_gamma(n_start) = [-0.05_dp, -0.10_dp, 0.00_dp, -0.15_dp]
        real(dp), parameter :: start_beta(n_start)  = [0.90_dp, 0.95_dp, 0.80_dp, 0.70_dp]
        real(dp) :: p(regarch1_np), p0(regarch1_np), p_best(regarch1_np), f_try
        real(dp) :: lv0, omega0
        integer :: istart, niter_try
        logical :: converged_try

        call regarch1_set_data(y, log_range)
        lv0 = regarch1_initial_log_volatility(log_range)
        f_best = huge(1.0_dp)
        p_best = 0.0_dp
        niter_best = 0
        converged_best = .false.

        do istart = 1, n_start
            omega0 = (1.0_dp - start_beta(istart)) * lv0
            call regarch1_inv_transform(omega0, start_alpha(istart), start_gamma(istart), start_beta(istart), p0)
            p = p0
            call bfgs_minimize(regarch1_obj, p, regarch1_np, max_iter, gtol, f_try, niter_try, converged_try)
            if (f_try < f_best) then
                f_best = f_try
                p_best = p
                niter_best = niter_try
                converged_best = converged_try
            end if
        end do

        params = garch_params_t()
        call regarch1_transform(p_best, params%omega, params%alpha, params%gamma, params%beta)
    end subroutine fit_regarch1

    subroutine fit_regarch2(y, log_range, max_iter, gtol, f_best, params, niter_best, converged_best)
        real(dp), intent(in)  :: y(:), log_range(:), gtol
        integer,  intent(in)  :: max_iter
        real(dp), intent(out) :: f_best
        type(garch_params_t), intent(out) :: params
        integer,  intent(out) :: niter_best
        logical,  intent(out) :: converged_best
        real(dp), parameter :: start_alpha_h(n_start) = [0.10_dp, 0.20_dp, 0.05_dp, 0.15_dp]
        real(dp), parameter :: start_gamma_h(n_start) = [-0.05_dp, -0.10_dp, 0.00_dp, -0.15_dp]
        real(dp), parameter :: start_beta_h(n_start)  = [0.80_dp, 0.90_dp, 0.70_dp, 0.85_dp]
        real(dp), parameter :: start_alpha_q(n_start) = [0.02_dp, 0.05_dp, 0.00_dp, 0.03_dp]
        real(dp), parameter :: start_gamma_q(n_start) = [0.00_dp, -0.02_dp, 0.00_dp, -0.05_dp]
        real(dp), parameter :: start_beta_q(n_start)  = [0.98_dp, 0.95_dp, 0.99_dp, 0.97_dp]
        real(dp) :: p(regarch2_np), p0(regarch2_np), p_best(regarch2_np), f_try
        real(dp) :: mu0
        integer :: istart, niter_try
        logical :: converged_try

        call regarch1_set_data(y, log_range)
        mu0 = regarch1_initial_log_volatility(log_range)
        f_best = huge(1.0_dp)
        p_best = 0.0_dp
        niter_best = 0
        converged_best = .false.

        do istart = 1, n_start
            call regarch2_inv_transform(mu0, start_alpha_h(istart), start_gamma_h(istart), &
                                        start_beta_h(istart), start_alpha_q(istart), &
                                        start_gamma_q(istart), start_beta_q(istart), p0)
            p = p0
            call bfgs_minimize(regarch2_obj, p, regarch2_np, max_iter, gtol, f_try, niter_try, converged_try)
            if (f_try < f_best) then
                f_best = f_try
                p_best = p
                niter_best = niter_try
                converged_best = converged_try
            end if
        end do

        params = garch_params_t()
        call regarch2_transform(p_best, params%omega, params%alpha, params%gamma, params%beta, &
                                params%theta, params%twist, params%scale)
    end subroutine fit_regarch2

    subroutine fit_rgarch_meas(y, log_range, max_iter, gtol, f_best, params, niter_best, converged_best)
        real(dp), intent(in)  :: y(:), log_range(:), gtol
        integer,  intent(in)  :: max_iter
        real(dp), intent(out) :: f_best
        type(garch_params_t), intent(out) :: params
        integer,  intent(out) :: niter_best
        logical,  intent(out) :: converged_best
        real(dp), parameter :: start_alpha(n_start) = [0.10_dp, 0.20_dp, 0.05_dp, 0.15_dp]
        real(dp), parameter :: start_gamma(n_start) = [-0.05_dp, -0.10_dp, 0.00_dp, -0.15_dp]
        real(dp), parameter :: start_beta(n_start)  = [0.90_dp, 0.95_dp, 0.80_dp, 0.98_dp]
        real(dp) :: p(rgarch_meas_np), p0(rgarch_meas_np), p_best(rgarch_meas_np), f_try
        real(dp) :: lv0, xi0
        integer :: istart, niter_try
        logical :: converged_try

        call regarch1_set_data(y, log_range)
        lv0 = log(max(sum(y**2) / real(size(y), dp), 1.0e-12_dp)) / 2.0_dp
        xi0 = sum(log_range) / real(size(log_range), dp) - lv0
        f_best = huge(1.0_dp)
        p_best = 0.0_dp
        niter_best = 0
        converged_best = .false.

        do istart = 1, n_start
            call rgarch_meas_inv_transform((1.0_dp - start_beta(istart))*lv0, start_alpha(istart), &
                                           start_gamma(istart), start_beta(istart), xi0, 1.0_dp, &
                                           0.0_dp, 0.0_dp, regarch_log_range_sd, p0)
            p = p0
            call bfgs_minimize(rgarch_meas_obj, p, rgarch_meas_np, max_iter, gtol, f_try, niter_try, converged_try)
            if (f_try < f_best) then
                f_best = f_try
                p_best = p
                niter_best = niter_try
                converged_best = converged_try
            end if
        end do

        params = garch_params_t()
        call rgarch_meas_transform(p_best, params%omega, params%alpha, params%gamma, params%beta, &
                                   params%theta, params%twist, params%extra1, params%extra2, params%scale)
    end subroutine fit_rgarch_meas

    subroutine fit_nagarch_range(y, range_var, max_iter, gtol, f_best, params, niter_best, converged_best)
        real(dp), intent(in)  :: y(:), range_var(:), gtol
        integer,  intent(in)  :: max_iter
        real(dp), intent(out) :: f_best
        type(garch_params_t), intent(out) :: params
        integer,  intent(out) :: niter_best
        logical,  intent(out) :: converged_best
        real(dp), parameter :: starts(nagarch_range_np,n_start) = reshape( &
            [1.0e-6_dp, 0.06_dp, 0.03_dp, 0.86_dp, 0.25_dp, &
             1.0e-6_dp, 0.08_dp, 0.05_dp, 0.82_dp, 0.50_dp, &
             5.0e-6_dp, 0.08_dp, 0.08_dp, 0.78_dp, 0.75_dp, &
             5.0e-6_dp, 0.10_dp, 0.05_dp, 0.78_dp, 1.00_dp], [nagarch_range_np,n_start])
        real(dp) :: p(nagarch_range_np), p0(nagarch_range_np), p_best(nagarch_range_np), f_try
        integer :: istart, niter_try
        logical :: converged_try

        call nagarch_range_set_data(y, range_var)
        f_best = huge(1.0_dp)
        p_best = 0.0_dp
        niter_best = 0
        converged_best = .false.

        do istart = 1, n_start
            call nagarch_range_inv_transform(starts(1,istart), starts(2,istart), starts(3,istart), &
                                             starts(4,istart), starts(5,istart), p0)
            p = p0
            call bfgs_minimize(nagarch_range_obj, p, nagarch_range_np, max_iter, gtol, f_try, niter_try, converged_try)
            if (f_try < f_best) then
                f_best = f_try
                p_best = p
                niter_best = niter_try
                converged_best = converged_try
            end if
        end do

        params = garch_params_t()
        call nagarch_range_transform(p_best, params%omega, params%alpha, params%gamma, params%beta, params%theta)
    end subroutine fit_nagarch_range

    subroutine fit_fgarch_twist_range(y, range_var, max_iter, gtol, f_best, params, niter_best, converged_best)
        real(dp), intent(in)  :: y(:), range_var(:), gtol
        integer,  intent(in)  :: max_iter
        real(dp), intent(out) :: f_best
        type(garch_params_t), intent(out) :: params
        integer,  intent(out) :: niter_best
        logical,  intent(out) :: converged_best
        real(dp), parameter :: starts(fgarch_twist_range_np,n_start) = reshape( &
            [1.0e-6_dp, 0.06_dp, 0.03_dp, 0.84_dp, 0.25_dp, 0.00_dp, &
             1.0e-6_dp, 0.08_dp, 0.05_dp, 0.80_dp, 0.50_dp, 0.00_dp, &
             5.0e-6_dp, 0.08_dp, 0.08_dp, 0.76_dp, 0.75_dp, 0.10_dp, &
             5.0e-6_dp, 0.10_dp, 0.05_dp, 0.76_dp, 1.00_dp, -0.10_dp], [fgarch_twist_range_np,n_start])
        real(dp) :: p(fgarch_twist_range_np), p0(fgarch_twist_range_np), p_best(fgarch_twist_range_np), f_try
        integer :: istart, niter_try
        logical :: converged_try

        call nagarch_range_set_data(y, range_var)
        f_best = huge(1.0_dp)
        p_best = 0.0_dp
        niter_best = 0
        converged_best = .false.

        do istart = 1, n_start
            call fgarch_twist_range_inv_transform(starts(1,istart), starts(2,istart), starts(3,istart), &
                                                  starts(4,istart), starts(5,istart), starts(6,istart), p0)
            p = p0
            call bfgs_minimize(fgarch_twist_range_obj, p, fgarch_twist_range_np, max_iter, gtol, &
                               f_try, niter_try, converged_try)
            if (f_try < f_best) then
                f_best = f_try
                p_best = p
                niter_best = niter_try
                converged_best = converged_try
            end if
        end do

        params = garch_params_t()
        call fgarch_twist_range_transform(p_best, params%omega, params%alpha, params%gamma, params%beta, &
                                          params%theta, params%twist)
    end subroutine fit_fgarch_twist_range

    subroutine fit_gjr(y, max_iter, gtol, f_best, params, niter_best, converged_best)
        real(dp), intent(in)  :: y(:), gtol
        integer,  intent(in)  :: max_iter
        real(dp), intent(out) :: f_best
        type(garch_params_t), intent(out) :: params
        integer,  intent(out) :: niter_best
        logical,  intent(out) :: converged_best
        real(dp), parameter :: starts(gjr_np,n_start) = reshape( &
            [1.0e-6_dp, 0.04_dp, 0.04_dp, 0.90_dp, &
             1.0e-6_dp, 0.04_dp, 0.10_dp, 0.88_dp, &
             1.0e-6_dp, 0.06_dp, 0.15_dp, 0.84_dp, &
             5.0e-6_dp, 0.08_dp, 0.08_dp, 0.82_dp], [gjr_np,n_start])
        real(dp) :: p(gjr_np), p0(gjr_np), p_best(gjr_np), f_try
        integer :: istart, niter_try, n
        logical :: converged_try

        n = size(y)
        call gjr_set_data(y, n)
        f_best = huge(1.0_dp)
        p_best = 0.0_dp
        niter_best = 0
        converged_best = .false.

        do istart = 1, n_start
            call gjr_inv_transform(starts(1,istart), starts(2,istart), &
                                   starts(3,istart), starts(4,istart), p0)
            p = p0
            call bfgs_minimize(gjr_obj, p, gjr_np, max_iter, gtol, f_try, niter_try, converged_try)
            if (f_try < f_best) then
                f_best = f_try
                p_best = p
                niter_best = niter_try
                converged_best = converged_try
            end if
        end do

        params = garch_params_t()
        call gjr_transform(p_best, params%omega, params%alpha, params%gamma, params%beta)
    end subroutine fit_gjr

    subroutine fit_gjr_signed(y, max_iter, gtol, f_best, params, niter_best, converged_best)
        real(dp), intent(in)  :: y(:), gtol
        integer,  intent(in)  :: max_iter
        real(dp), intent(out) :: f_best
        type(garch_params_t), intent(out) :: params
        integer,  intent(out) :: niter_best
        logical,  intent(out) :: converged_best
        real(dp), parameter :: starts(gjr_np,n_start) = reshape( &
            [1.0e-6_dp, 0.04_dp,  0.04_dp, 0.90_dp, &
             1.0e-6_dp, 0.04_dp,  0.10_dp, 0.88_dp, &
             1.0e-6_dp, 0.06_dp, -0.02_dp, 0.93_dp, &
             5.0e-6_dp, 0.08_dp, -0.04_dp, 0.88_dp], [gjr_np,n_start])
        real(dp) :: p(gjr_np), p0(gjr_np), p_best(gjr_np), f_try
        integer :: istart, niter_try, n
        logical :: converged_try

        n = size(y)
        call gjr_set_data(y, n)
        f_best = huge(1.0_dp)
        p_best = 0.0_dp
        niter_best = 0
        converged_best = .false.

        do istart = 1, n_start
            call gjr_signed_inv_transform(starts(1,istart), starts(2,istart), &
                                          starts(3,istart), starts(4,istart), p0)
            p = p0
            call bfgs_minimize(gjr_signed_obj, p, gjr_np, max_iter, gtol, f_try, niter_try, converged_try)
            if (f_try < f_best) then
                f_best = f_try
                p_best = p
                niter_best = niter_try
                converged_best = converged_try
            end if
        end do

        params = garch_params_t()
        call gjr_signed_transform(p_best, params%omega, params%alpha, params%gamma, params%beta)
    end subroutine fit_gjr_signed

    subroutine fit_egarch(y, max_iter, gtol, f_best, params, niter_best, converged_best)
        real(dp), intent(in)  :: y(:), gtol
        integer,  intent(in)  :: max_iter
        real(dp), intent(out) :: f_best
        type(garch_params_t), intent(out) :: params
        integer,  intent(out) :: niter_best
        logical,  intent(out) :: converged_best
        real(dp), parameter :: start_alpha(n_start) = [0.08_dp, 0.10_dp, 0.06_dp, 0.12_dp]
        real(dp), parameter :: start_gamma(n_start) = [-0.05_dp, -0.10_dp, 0.00_dp, -0.15_dp]
        real(dp), parameter :: start_beta(n_start)  = [0.90_dp, 0.95_dp, 0.85_dp, 0.80_dp]
        real(dp) :: p(egarch_np), p0(egarch_np), p_best(egarch_np), f_try
        real(dp) :: target_lh, omega0
        integer :: istart, niter_try, n
        logical :: converged_try

        n = size(y)
        call egarch_set_data(y, n)
        target_lh = log(max(sum(y**2) / real(n, dp), 1.0e-12_dp))
        f_best = huge(1.0_dp)
        p_best = 0.0_dp
        niter_best = 0
        converged_best = .false.

        do istart = 1, n_start
            omega0 = target_lh * (1.0_dp - start_beta(istart))
            call egarch_inv_transform(omega0, start_alpha(istart), start_gamma(istart), start_beta(istart), p0)
            p = p0
            call bfgs_minimize(egarch_obj, p, egarch_np, max_iter, gtol, f_try, niter_try, converged_try)
            if (f_try < f_best) then
                f_best = f_try
                p_best = p
                niter_best = niter_try
                converged_best = converged_try
            end if
        end do

        params = garch_params_t()
        call egarch_transform(p_best, params%omega, params%alpha, params%gamma, params%beta)
    end subroutine fit_egarch

    subroutine fit_aparch(y, max_iter, gtol, f_best, params, niter_best, converged_best)
        real(dp), intent(in)  :: y(:), gtol
        integer,  intent(in)  :: max_iter
        real(dp), intent(out) :: f_best
        type(garch_params_t), intent(out) :: params
        integer,  intent(out) :: niter_best
        logical,  intent(out) :: converged_best
        real(dp), parameter :: starts(aparch_np,n_start) = reshape( &
            [1.0e-4_dp, 0.08_dp,  0.50_dp, 0.88_dp, 1.00_dp, &
             1.0e-4_dp, 0.10_dp,  0.80_dp, 0.86_dp, 1.00_dp, &
             1.0e-4_dp, 0.08_dp, -0.20_dp, 0.88_dp, 1.50_dp, &
             5.0e-5_dp, 0.05_dp,  0.50_dp, 0.92_dp, 2.00_dp], [aparch_np,n_start])
        real(dp) :: p(aparch_np), p0(aparch_np), p_best(aparch_np), f_try
        real(dp) :: omega0
        integer :: istart, niter_try, n
        logical :: converged_try

        n = size(y)
        call aparch_set_data(y)
        f_best = huge(1.0_dp)
        p_best = 0.0_dp
        niter_best = 0
        converged_best = .false.

        do istart = 1, n_start
            omega0 = max((1.0_dp - starts(2,istart) - starts(4,istart)) * &
                         sum(abs(y)**starts(5,istart)) / real(n, dp), 1.0e-12_dp)
            call aparch_inv_transform(omega0, starts(2,istart), starts(3,istart), &
                                      starts(4,istart), starts(5,istart), p0)
            p = p0
            call bfgs_minimize(aparch_obj, p, aparch_np, max_iter, gtol, f_try, niter_try, converged_try)
            if (f_try < f_best) then
                f_best = f_try
                p_best = p
                niter_best = niter_try
                converged_best = converged_try
            end if
        end do

        params = garch_params_t()
        call aparch_transform(p_best, params%omega, params%alpha, params%gamma, params%beta, params%theta)
    end subroutine fit_aparch

    subroutine fit_fgarch_twist(y, ret_std, max_iter, gtol, f_best, params, vol_ann_best, skew_best, ekurt_best, &
                                niter_best, converged_best)
        real(dp), intent(in)  :: y(:), ret_std, gtol
        integer,  intent(in)  :: max_iter
        real(dp), intent(out) :: f_best
        type(garch_params_t), intent(out) :: params
        real(dp), intent(out) :: vol_ann_best, skew_best, ekurt_best
        integer,  intent(out) :: niter_best
        logical,  intent(out) :: converged_best
        integer :: n
        real(dp) :: omega0
        real(dp) :: p(fgarch_twist_np), p0(fgarch_twist_np), p_best(fgarch_twist_np)

        call fgarch_set_dist(fg_dist_normal)
        n = size(y)
        call fgarch_set_data(y, n)
        omega0 = max(1.0e-8_dp, ret_std**2 * (1.0_dp - 0.90_dp - 0.05_dp))
        call fgarch_twist_inv_transform(omega0, 0.05_dp, 0.90_dp, 0.0_dp, 0.0_dp, p0)
        p = p0
        call bfgs_minimize(fgarch_twist_obj, p, fgarch_twist_np, max_iter, gtol, f_best, niter_best, converged_best)
        p_best = p

        params = garch_params_t()
        call fgarch_twist_transform(p_best, params%omega, params%alpha, params%beta, params%theta, params%twist)
        call fgarch_twist_vol_ann(p_best, fgarch_twist_np, vol_ann_best)
        call fgarch_twist_skew_kurt(p_best, fgarch_twist_np, skew_best, ekurt_best)
    end subroutine fit_fgarch_twist

    subroutine aparch_set_data(y)
        real(dp), intent(in) :: y(:)
        if (allocated(aparch_obs)) deallocate(aparch_obs)
        allocate(aparch_obs(size(y)))
        aparch_obs = y
        aparch_nobs = size(y)
    end subroutine aparch_set_data

    subroutine harch_set_data(y)
        real(dp), intent(in) :: y(:)
        if (allocated(harch_obs)) deallocate(harch_obs)
        allocate(harch_obs(size(y)))
        harch_obs = y
        harch_nobs = size(y)
    end subroutine harch_set_data

    subroutine qgarch_set_data(y)
        real(dp), intent(in) :: y(:)
        if (allocated(qgarch_obs)) deallocate(qgarch_obs)
        allocate(qgarch_obs(size(y)))
        qgarch_obs = y
        qgarch_nobs = size(y)
    end subroutine qgarch_set_data

    subroutine figarch_set_data(y)
        real(dp), intent(in) :: y(:)
        if (allocated(figarch_obs)) deallocate(figarch_obs)
        allocate(figarch_obs(size(y)))
        figarch_obs = y
        figarch_nobs = size(y)
    end subroutine figarch_set_data

    subroutine qgarch_transform(p, omega, alpha, beta, shift)
        real(dp), intent(in)  :: p(qgarch_np)
        real(dp), intent(out) :: omega, alpha, beta, shift
        real(dp) :: ea, eb, s

        omega = exp(p(1))
        ea = exp(p(2))
        eb = exp(p(3))
        s = 1.0_dp + ea + eb
        alpha = ea / s
        beta = eb / s
        shift = p(4)
    end subroutine qgarch_transform

    subroutine qgarch_inv_transform(omega, alpha, beta, shift, p)
        real(dp), intent(in)  :: omega, alpha, beta, shift
        real(dp), intent(out) :: p(qgarch_np)
        real(dp) :: slack

        slack = max(1.0_dp - alpha - beta, 1.0e-8_dp)
        p(1) = log(max(omega, 1.0e-12_dp))
        p(2) = log(max(alpha, 1.0e-12_dp) / slack)
        p(3) = log(max(beta, 1.0e-12_dp) / slack)
        p(4) = shift
    end subroutine qgarch_inv_transform

    real(dp) function qgarch_initial_variance(omega, alpha, beta, shift)
        real(dp), intent(in) :: omega, alpha, beta, shift

        qgarch_initial_variance = max((omega + alpha*shift**2) / max(1.0_dp - alpha - beta, 1.0e-8_dp), &
                                      1.0e-12_dp)
    end function qgarch_initial_variance

    subroutine qgarch_obj(p, np, f, g)
        integer,  intent(in)  :: np
        real(dp), intent(in)  :: p(np)
        real(dp), intent(out) :: f
        real(dp), intent(out) :: g(np)
        real(dp) :: omega, alpha, beta, shift, h, h_old, y, x, factor
        real(dp) :: den, numer
        real(dp) :: dh_om, dh_al, dh_be, dh_sh
        real(dp) :: dh_om_old, dh_al_old, dh_be_old, dh_sh_old
        real(dp) :: grad_om, grad_al, grad_be, grad_sh, weighted
        integer :: t

        call qgarch_transform(p, omega, alpha, beta, shift)
        den = max(1.0_dp - alpha - beta, 1.0e-8_dp)
        numer = omega + alpha*shift**2
        h = max(numer / den, 1.0e-12_dp)
        dh_om = 1.0_dp / den
        dh_al = (shift**2*den + numer) / den**2
        dh_be = numer / den**2
        dh_sh = 2.0_dp*alpha*shift / den

        f = real(qgarch_nobs, dp) * log_sqrt_2pi
        grad_om = 0.0_dp
        grad_al = 0.0_dp
        grad_be = 0.0_dp
        grad_sh = 0.0_dp
        do t = 1, qgarch_nobs
            h = max(h, 1.0e-12_dp)
            y = qgarch_obs(t)
            factor = 1.0_dp/h - y**2 / h**2
            f = f + 0.5_dp * (log(h) + y**2 / h)
            grad_om = grad_om + 0.5_dp * factor * dh_om
            grad_al = grad_al + 0.5_dp * factor * dh_al
            grad_be = grad_be + 0.5_dp * factor * dh_be
            grad_sh = grad_sh + 0.5_dp * factor * dh_sh

            h_old = h
            dh_om_old = dh_om
            dh_al_old = dh_al
            dh_be_old = dh_be
            dh_sh_old = dh_sh
            x = y - shift
            h = omega + alpha*x**2 + beta*h_old
            dh_om = 1.0_dp + beta*dh_om_old
            dh_al = x**2 + beta*dh_al_old
            dh_be = h_old + beta*dh_be_old
            dh_sh = -2.0_dp*alpha*x + beta*dh_sh_old
        end do

        weighted = alpha*grad_al + beta*grad_be
        g(1) = grad_om * omega
        g(2) = alpha * (grad_al - weighted)
        g(3) = beta * (grad_be - weighted)
        g(4) = grad_sh

        f = f / real(qgarch_nobs, dp)
        g = g / real(qgarch_nobs, dp)
    end subroutine qgarch_obj

    logical function qgarch_nagarch_start(y, max_iter, gtol, p0)
        real(dp), intent(in)  :: y(:), gtol
        integer,  intent(in)  :: max_iter
        real(dp), intent(out) :: p0(qgarch_np)
        type(garch_params_t) :: nag_params
        real(dp) :: f_nag, h, sqrth, r, mean_h, mean_sqrth, q_shift
        real(dp) :: alpha0, beta0, omega0, slack, var_y
        integer :: t, niter_nag
        logical :: conv_nag

        qgarch_nagarch_start = .false.
        p0 = 0.0_dp

        call fit_nagarch(y, max_iter, gtol, f_nag, nag_params, niter_nag, conv_nag)
        if (.not. conv_nag) return

        h = nag_params%omega / max(1.0_dp - nagarch_persist(nag_params), 1.0e-8_dp)
        h = max(h, 1.0e-12_dp)
        mean_h = 0.0_dp
        mean_sqrth = 0.0_dp
        do t = 1, size(y)
            h = max(h, 1.0e-12_dp)
            sqrth = sqrt(h)
            mean_h = mean_h + h
            mean_sqrth = mean_sqrth + sqrth
            r = y(t) - nag_params%theta*sqrth
            h = nag_params%omega + nag_params%alpha*r**2 + nag_params%beta*h
        end do
        mean_h = mean_h / real(size(y), dp)
        mean_sqrth = mean_sqrth / real(size(y), dp)

        alpha0 = nag_params%alpha
        beta0 = nag_params%beta
        slack = 1.0_dp - alpha0 - beta0
        if (slack <= 1.0e-6_dp) then
            alpha0 = min(alpha0, 0.95_dp*(1.0_dp - beta0))
            slack = 1.0_dp - alpha0 - beta0
        end if
        if (slack <= 1.0e-6_dp) return

        q_shift = nag_params%theta * mean_sqrth
        omega0 = mean_h*slack - alpha0*q_shift**2
        if (omega0 <= 1.0e-12_dp) then
            var_y = max(sum(y**2) / real(size(y), dp), 1.0e-12_dp)
            omega0 = max(0.10_dp*slack*var_y, 1.0e-12_dp)
        end if

        call qgarch_inv_transform(omega0, alpha0, beta0, q_shift, p0)
        qgarch_nagarch_start = .true.
    end function qgarch_nagarch_start

    subroutine figarch_transform(p, omega, phi, d, beta)
        real(dp), intent(in)  :: p(figarch_np)
        real(dp), intent(out) :: omega, phi, d, beta
        real(dp) :: ud, uphi, ubeta, phi_max

        omega = exp(p(1))
        ud = 1.0_dp / (1.0_dp + exp(-p(2)))
        uphi = 1.0_dp / (1.0_dp + exp(-p(3)))
        d = ud
        phi_max = 0.5_dp * (1.0_dp - d)
        phi = phi_max * uphi
        ubeta = 1.0_dp / (1.0_dp + exp(-p(4)))
        beta = (d + phi) * ubeta
    end subroutine figarch_transform

    subroutine figarch_inv_transform(omega, phi, d, beta, p)
        real(dp), intent(in)  :: omega, phi, d, beta
        real(dp), intent(out) :: p(figarch_np)
        real(dp) :: dd, pp, bb, phi_max

        dd = min(max(d, 1.0e-8_dp), 1.0_dp - 1.0e-8_dp)
        phi_max = max(0.5_dp * (1.0_dp - dd), 1.0e-8_dp)
        pp = min(max(phi / phi_max, 1.0e-8_dp), 1.0_dp - 1.0e-8_dp)
        bb = min(max(beta / max(dd + phi, 1.0e-8_dp), 1.0e-8_dp), 1.0_dp - 1.0e-8_dp)
        p(1) = log(max(omega, 1.0e-12_dp))
        p(2) = log(dd / (1.0_dp - dd))
        p(3) = log(pp / (1.0_dp - pp))
        p(4) = log(bb / (1.0_dp - bb))
    end subroutine figarch_inv_transform

    subroutine fi_nagarch_transform(p, omega, phi, d, beta, shift)
        real(dp), intent(in)  :: p(fi_nagarch_np)
        real(dp), intent(out) :: omega, phi, d, beta, shift

        call figarch_transform(p(1:figarch_np), omega, phi, d, beta)
        shift = p(5)
    end subroutine fi_nagarch_transform

    subroutine fi_nagarch_inv_transform(omega, phi, d, beta, shift, p)
        real(dp), intent(in)  :: omega, phi, d, beta, shift
        real(dp), intent(out) :: p(fi_nagarch_np)

        call figarch_inv_transform(omega, phi, d, beta, p(1:figarch_np))
        p(5) = shift
    end subroutine fi_nagarch_inv_transform

    subroutine figarch_weights(phi, d, beta, lambda)
        real(dp), intent(in)  :: phi, d, beta
        real(dp), intent(out) :: lambda(:)
        real(dp) :: delta_prev, delta_cur
        integer :: i

        if (size(lambda) < 1) return
        lambda(1) = phi - beta + d
        delta_prev = d
        do i = 2, size(lambda)
            delta_cur = (real(i - 1, dp) - d) / real(i, dp) * delta_prev
            lambda(i) = beta*lambda(i - 1) + delta_cur - phi*delta_prev
            delta_prev = delta_cur
        end do
    end subroutine figarch_weights

    subroutine figarch_weights_deriv(phi, d, beta, lambda, dl_dphi, dl_dd, dl_dbeta)
        real(dp), intent(in)  :: phi, d, beta
        real(dp), intent(out) :: lambda(:), dl_dphi(:), dl_dd(:), dl_dbeta(:)
        real(dp) :: delta_prev, delta_cur, ddelta_prev, ddelta_cur, coeff
        integer :: i

        if (size(lambda) < 1) return
        lambda(1) = phi - beta + d
        dl_dphi(1) = 1.0_dp
        dl_dd(1) = 1.0_dp
        dl_dbeta(1) = -1.0_dp
        delta_prev = d
        ddelta_prev = 1.0_dp
        do i = 2, size(lambda)
            coeff = (real(i - 1, dp) - d) / real(i, dp)
            delta_cur = coeff * delta_prev
            ddelta_cur = -delta_prev / real(i, dp) + coeff * ddelta_prev
            lambda(i) = beta*lambda(i - 1) + delta_cur - phi*delta_prev
            dl_dphi(i) = beta*dl_dphi(i - 1) - delta_prev
            dl_dd(i) = beta*dl_dd(i - 1) + ddelta_cur - phi*ddelta_prev
            dl_dbeta(i) = lambda(i - 1) + beta*dl_dbeta(i - 1)
            delta_prev = delta_cur
            ddelta_prev = ddelta_cur
        end do
    end subroutine figarch_weights_deriv

    real(dp) function figarch_backcast(y)
        real(dp), intent(in) :: y(:)
        real(dp) :: w, wsum
        integer :: i, tau

        tau = min(75, size(y))
        figarch_backcast = 0.0_dp
        wsum = 0.0_dp
        do i = 1, tau
            w = 0.94_dp**real(i - 1, dp)
            figarch_backcast = figarch_backcast + w*y(i)**2
            wsum = wsum + w
        end do
        figarch_backcast = max(figarch_backcast / max(wsum, 1.0e-12_dp), 1.0e-12_dp)
    end function figarch_backcast

    subroutine figarch_variance(y, params, variance)
        real(dp), intent(in) :: y(:)
        type(garch_params_t), intent(in) :: params
        real(dp), intent(out) :: variance(:)
        real(dp), allocatable :: lambda(:)
        real(dp) :: omega_tilde, backcast, bc_weight
        integer :: t, i, n, m

        n = size(y)
        m = figarch_trunc_lag
        allocate(lambda(m))
        call figarch_weights(params%alpha, params%theta, params%beta, lambda)
        omega_tilde = params%omega / max(1.0_dp - params%beta, 1.0e-8_dp)
        backcast = figarch_backcast(y)

        do t = 1, n
            bc_weight = 0.0_dp
            do i = t, m
                bc_weight = bc_weight + lambda(i)
            end do
            variance(t) = omega_tilde + bc_weight*backcast
            do i = 1, min(t - 1, m)
                variance(t) = variance(t) + lambda(i)*y(t - i)**2
            end do
            variance(t) = max(variance(t), 1.0e-12_dp)
        end do
        deallocate(lambda)
    end subroutine figarch_variance

    subroutine fi_nagarch_variance(y, params, variance)
        real(dp), intent(in) :: y(:)
        type(garch_params_t), intent(in) :: params
        real(dp), intent(out) :: variance(:)
        real(dp), allocatable :: lambda(:), news(:)
        real(dp) :: omega_tilde, backcast, bc_weight, h, sqrth, scale
        integer :: t, i, n, m

        n = size(y)
        m = figarch_trunc_lag
        allocate(lambda(m), news(n))
        call figarch_weights(params%alpha, params%theta, params%beta, lambda)
        omega_tilde = params%omega / max(1.0_dp - params%beta, 1.0e-8_dp)
        backcast = figarch_backcast(y)
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
            variance(t) = h
            sqrth = sqrt(h)
            news(t) = (y(t) - params%twist*sqrth)**2 / scale
        end do
        deallocate(lambda, news)
    end subroutine fi_nagarch_variance

    real(dp) function figarch_value(p)
        real(dp), intent(in) :: p(figarch_np)
        type(garch_params_t) :: params
        real(dp), allocatable :: variance(:)
        integer :: t

        call figarch_transform(p, params%omega, params%alpha, params%theta, params%beta)
        allocate(variance(figarch_nobs))
        call figarch_variance(figarch_obs, params, variance)
        figarch_value = real(figarch_nobs, dp) * log_sqrt_2pi
        do t = 1, figarch_nobs
            figarch_value = figarch_value + 0.5_dp * (log(variance(t)) + figarch_obs(t)**2 / variance(t))
        end do
        figarch_value = figarch_value / real(figarch_nobs, dp)
        deallocate(variance)
    end function figarch_value

    real(dp) function fi_nagarch_value(p)
        real(dp), intent(in) :: p(fi_nagarch_np)
        type(garch_params_t) :: params
        real(dp), allocatable :: variance(:)
        integer :: t

        call fi_nagarch_transform(p, params%omega, params%alpha, params%theta, params%beta, params%twist)
        allocate(variance(figarch_nobs))
        call fi_nagarch_variance(figarch_obs, params, variance)
        fi_nagarch_value = real(figarch_nobs, dp) * log_sqrt_2pi
        do t = 1, figarch_nobs
            fi_nagarch_value = fi_nagarch_value + 0.5_dp * (log(variance(t)) + figarch_obs(t)**2 / variance(t))
        end do
        fi_nagarch_value = fi_nagarch_value / real(figarch_nobs, dp)
        deallocate(variance)
    end function fi_nagarch_value

    subroutine figarch_obj(p, np, f, g)
        integer,  intent(in)  :: np
        real(dp), intent(in)  :: p(np)
        real(dp), intent(out) :: f
        real(dp), intent(out) :: g(np)
        real(dp), allocatable :: lambda(:), dl_dphi(:), dl_dd(:), dl_dbeta(:)
        real(dp) :: omega, phi, d, beta, omega_tilde, backcast, bc_weight
        real(dp) :: h, y, factor, lag_sq
        real(dp) :: dh_om, dh_phi, dh_d, dh_beta
        real(dp) :: grad_om, grad_phi, grad_d, grad_beta
        real(dp) :: ud, uphi, ubeta, phi_max
        real(dp) :: dd_dp2, dphi_dp2, dphi_dp3, dbeta_dp2, dbeta_dp3, dbeta_dp4
        integer :: t, i, m

        call figarch_transform(p, omega, phi, d, beta)
        m = figarch_trunc_lag
        allocate(lambda(m), dl_dphi(m), dl_dd(m), dl_dbeta(m))
        call figarch_weights_deriv(phi, d, beta, lambda, dl_dphi, dl_dd, dl_dbeta)
        omega_tilde = omega / max(1.0_dp - beta, 1.0e-8_dp)
        backcast = figarch_backcast(figarch_obs)

        f = real(figarch_nobs, dp) * log_sqrt_2pi
        grad_om = 0.0_dp
        grad_phi = 0.0_dp
        grad_d = 0.0_dp
        grad_beta = 0.0_dp
        do t = 1, figarch_nobs
            bc_weight = 0.0_dp
            dh_phi = 0.0_dp
            dh_d = 0.0_dp
            dh_beta = omega / max(1.0_dp - beta, 1.0e-8_dp)**2
            do i = t, m
                bc_weight = bc_weight + lambda(i)
                dh_phi = dh_phi + dl_dphi(i)*backcast
                dh_d = dh_d + dl_dd(i)*backcast
                dh_beta = dh_beta + dl_dbeta(i)*backcast
            end do
            h = omega_tilde + bc_weight*backcast
            do i = 1, min(t - 1, m)
                lag_sq = figarch_obs(t - i)**2
                h = h + lambda(i)*lag_sq
                dh_phi = dh_phi + dl_dphi(i)*lag_sq
                dh_d = dh_d + dl_dd(i)*lag_sq
                dh_beta = dh_beta + dl_dbeta(i)*lag_sq
            end do
            h = max(h, 1.0e-12_dp)
            y = figarch_obs(t)
            factor = 1.0_dp/h - y**2 / h**2
            f = f + 0.5_dp * (log(h) + y**2 / h)
            dh_om = 1.0_dp / max(1.0_dp - beta, 1.0e-8_dp)
            grad_om = grad_om + 0.5_dp * factor * dh_om
            grad_phi = grad_phi + 0.5_dp * factor * dh_phi
            grad_d = grad_d + 0.5_dp * factor * dh_d
            grad_beta = grad_beta + 0.5_dp * factor * dh_beta
        end do

        ud = d
        phi_max = max(0.5_dp * (1.0_dp - d), 1.0e-8_dp)
        uphi = min(max(phi / phi_max, 1.0e-8_dp), 1.0_dp - 1.0e-8_dp)
        ubeta = min(max(beta / max(d + phi, 1.0e-8_dp), 1.0e-8_dp), 1.0_dp - 1.0e-8_dp)
        dd_dp2 = ud * (1.0_dp - ud)
        dphi_dp2 = -0.5_dp * uphi * dd_dp2
        dphi_dp3 = phi_max * uphi * (1.0_dp - uphi)
        dbeta_dp2 = ubeta * (dd_dp2 + dphi_dp2)
        dbeta_dp3 = ubeta * dphi_dp3
        dbeta_dp4 = (d + phi) * ubeta * (1.0_dp - ubeta)

        g(1) = grad_om * omega
        g(2) = grad_d*dd_dp2 + grad_phi*dphi_dp2 + grad_beta*dbeta_dp2
        g(3) = grad_phi*dphi_dp3 + grad_beta*dbeta_dp3
        g(4) = grad_beta*dbeta_dp4

        f = f / real(figarch_nobs, dp)
        g = g / real(figarch_nobs, dp)
        deallocate(lambda, dl_dphi, dl_dd, dl_dbeta)
    end subroutine figarch_obj

    subroutine fi_nagarch_obj(p, np, f, g)
        integer,  intent(in)  :: np
        real(dp), intent(in)  :: p(np)
        real(dp), intent(out) :: f
        real(dp), intent(out) :: g(np)
        real(dp) :: pp(fi_nagarch_np), fp, fm, step
        integer :: i

        f = fi_nagarch_value(p)
        do i = 1, np
            step = 1.0e-5_dp * max(abs(p(i)), 1.0_dp)
            pp = p
            pp(i) = p(i) + step
            fp = fi_nagarch_value(pp)
            pp = p
            pp(i) = p(i) - step
            fm = fi_nagarch_value(pp)
            g(i) = (fp - fm) / (2.0_dp*step)
        end do
    end subroutine fi_nagarch_obj

    subroutine harch_transform(p, omega, alpha1, alpha5, alpha22)
        real(dp), intent(in)  :: p(harch_np)
        real(dp), intent(out) :: omega, alpha1, alpha5, alpha22
        real(dp) :: e2, e3, e4, s
        omega = exp(p(1))
        e2 = exp(p(2))
        e3 = exp(p(3))
        e4 = exp(p(4))
        s = 1.0_dp + e2 + e3 + e4
        alpha1 = e2 / s
        alpha5 = e3 / s
        alpha22 = e4 / s
    end subroutine harch_transform

    subroutine harch_inv_transform(omega, alpha1, alpha5, alpha22, p)
        real(dp), intent(in)  :: omega, alpha1, alpha5, alpha22
        real(dp), intent(out) :: p(harch_np)
        real(dp) :: slack
        slack = max(1.0_dp - alpha1 - alpha5 - alpha22, 1.0e-8_dp)
        p(1) = log(max(omega, 1.0e-12_dp))
        p(2) = log(max(alpha1, 1.0e-12_dp) / slack)
        p(3) = log(max(alpha5, 1.0e-12_dp) / slack)
        p(4) = log(max(alpha22, 1.0e-12_dp) / slack)
    end subroutine harch_inv_transform

    real(dp) function harch_block(t, lag, backcast)
        integer, intent(in) :: t, lag
        real(dp), intent(in) :: backcast
        integer :: j, idx
        harch_block = 0.0_dp
        do j = 1, lag
            idx = t - j
            if (idx >= 1) then
                harch_block = harch_block + harch_obs(idx)**2
            else
                harch_block = harch_block + backcast
            end if
        end do
        harch_block = harch_block / real(lag, dp)
    end function harch_block

    subroutine harch_obj(p, np, f, g)
        integer,  intent(in)  :: np
        real(dp), intent(in)  :: p(np)
        real(dp), intent(out) :: f
        real(dp), intent(out) :: g(np)
        real(dp) :: omega, alpha1, alpha5, alpha22, x1, x5, x22, h, y, factor, backcast
        real(dp) :: grad_om, grad_a1, grad_a5, grad_a22, weighted
        integer :: t

        call harch_transform(p, omega, alpha1, alpha5, alpha22)
        backcast = max(sum(harch_obs**2) / real(harch_nobs, dp), 1.0e-12_dp)

        f = real(harch_nobs, dp) * log_sqrt_2pi
        grad_om = 0.0_dp
        grad_a1 = 0.0_dp
        grad_a5 = 0.0_dp
        grad_a22 = 0.0_dp
        do t = 1, harch_nobs
            x1 = harch_block(t, 1, backcast)
            x5 = harch_block(t, 5, backcast)
            x22 = harch_block(t, 22, backcast)
            h = max(omega + alpha1*x1 + alpha5*x5 + alpha22*x22, 1.0e-12_dp)
            y = harch_obs(t)
            factor = 1.0_dp/h - y**2 / h**2
            f = f + 0.5_dp * (log(h) + y**2 / h)
            grad_om = grad_om + 0.5_dp * factor
            grad_a1 = grad_a1 + 0.5_dp * factor * x1
            grad_a5 = grad_a5 + 0.5_dp * factor * x5
            grad_a22 = grad_a22 + 0.5_dp * factor * x22
        end do

        weighted = alpha1*grad_a1 + alpha5*grad_a5 + alpha22*grad_a22
        g(1) = grad_om * omega
        g(2) = alpha1 * (grad_a1 - weighted)
        g(3) = alpha5 * (grad_a5 - weighted)
        g(4) = alpha22 * (grad_a22 - weighted)

        f = f / real(harch_nobs, dp)
        g = g / real(harch_nobs, dp)
    end subroutine harch_obj

    subroutine riskmetrics2006_settings(weights, mus)
        real(dp), intent(out) :: weights(rm2006_kmax), mus(rm2006_kmax)
        real(dp), parameter :: tau0 = 1560.0_dp, tau1 = 4.0_dp, rho = sqrt2
        real(dp) :: tau, wsum
        integer :: k

        wsum = 0.0_dp
        do k = 1, rm2006_kmax
            tau = tau1 * rho**real(k - 1, dp)
            weights(k) = 1.0_dp - log(tau) / log(tau0)
            mus(k) = exp(-1.0_dp / tau)
            wsum = wsum + weights(k)
        end do
        weights = weights / wsum
    end subroutine riskmetrics2006_settings

    subroutine riskmetrics2006_backcast(y, mus, backcast)
        real(dp), intent(in) :: y(:), mus(rm2006_kmax)
        real(dp), intent(out) :: backcast(rm2006_kmax)
        real(dp) :: weight_sum, w
        integer :: k, j, endpoint, n

        n = size(y)
        do k = 1, rm2006_kmax
            endpoint = int(max(min(real(floor(log(0.01_dp) / log(mus(k))), dp), real(n, dp)), real(k - 1, dp)))
            endpoint = max(endpoint, 1)
            weight_sum = 0.0_dp
            backcast(k) = 0.0_dp
            do j = 1, endpoint
                w = mus(k)**real(j - 1, dp)
                weight_sum = weight_sum + w
                backcast(k) = backcast(k) + w*y(j)**2
            end do
            backcast(k) = backcast(k) / weight_sum
        end do
    end subroutine riskmetrics2006_backcast

    subroutine riskmetrics2006_variance(y, variance)
        real(dp), intent(in) :: y(:)
        real(dp), intent(out) :: variance(:)
        real(dp) :: weights(rm2006_kmax), mus(rm2006_kmax), backcast(rm2006_kmax)
        real(dp) :: comp(rm2006_kmax)
        integer :: t, k, n

        n = size(y)
        if (size(variance) /= n) then
            print '(A)', "garch_fit_mod: variance length mismatch in riskmetrics2006_variance"
            error stop
        end if

        call riskmetrics2006_settings(weights, mus)
        call riskmetrics2006_backcast(y, mus, backcast)
        comp = backcast
        do t = 1, n
            variance(t) = sum(weights * comp)
            do k = 1, rm2006_kmax
                comp(k) = mus(k)*comp(k) + (1.0_dp - mus(k))*y(t)**2
            end do
        end do
    end subroutine riskmetrics2006_variance

    subroutine midas_hyperbolic_set_data(y)
        real(dp), intent(in) :: y(:)
        if (allocated(midas_hyperbolic_obs)) deallocate(midas_hyperbolic_obs)
        allocate(midas_hyperbolic_obs(size(y)))
        midas_hyperbolic_obs = y
        midas_hyperbolic_nobs = size(y)
    end subroutine midas_hyperbolic_set_data

    subroutine midas_hyperbolic_transform(p, omega, alpha, theta)
        real(dp), intent(in)  :: p(midas_hyperbolic_np)
        real(dp), intent(out) :: omega, alpha, theta
        omega = exp(p(1))
        alpha = 1.0_dp / (1.0_dp + exp(-p(2)))
        theta = 1.0_dp / (1.0_dp + exp(-p(3)))
    end subroutine midas_hyperbolic_transform

    subroutine midas_hyperbolic_inv_transform(omega, alpha, theta, p)
        real(dp), intent(in)  :: omega, alpha, theta
        real(dp), intent(out) :: p(midas_hyperbolic_np)
        real(dp) :: a, th
        a = min(max(alpha, 1.0e-8_dp), 1.0_dp - 1.0e-8_dp)
        th = min(max(theta, 1.0e-8_dp), 1.0_dp - 1.0e-8_dp)
        p(1) = log(max(omega, 1.0e-12_dp))
        p(2) = log(a / (1.0_dp - a))
        p(3) = log(th / (1.0_dp - th))
    end subroutine midas_hyperbolic_inv_transform

    subroutine midas_hyperbolic_asym_transform(p, omega, alpha, gamma, theta)
        real(dp), intent(in)  :: p(midas_hyperbolic_asym_np)
        real(dp), intent(out) :: omega, alpha, gamma, theta
        real(dp) :: amax

        omega = exp(p(1))
        gamma = -1.0_dp + 3.0_dp / (1.0_dp + exp(-p(3)))
        amax = min(1.0_dp - 0.5_dp*gamma, 1.0_dp)
        amax = max(amax, 1.0e-8_dp)
        alpha = amax / (1.0_dp + exp(-p(2)))
        theta = 1.0_dp / (1.0_dp + exp(-p(4)))
    end subroutine midas_hyperbolic_asym_transform

    subroutine midas_hyperbolic_asym_inv_transform(omega, alpha, gamma, theta, p)
        real(dp), intent(in)  :: omega, alpha, gamma, theta
        real(dp), intent(out) :: p(midas_hyperbolic_asym_np)
        real(dp) :: a, g, th, amax

        g = min(max(gamma, -1.0_dp + 1.0e-8_dp), 2.0_dp - 1.0e-8_dp)
        amax = min(1.0_dp - 0.5_dp*g, 1.0_dp)
        amax = max(amax, 1.0e-8_dp)
        a = min(max(alpha, max(-g + 1.0e-8_dp, 1.0e-8_dp)), amax - 1.0e-8_dp)
        th = min(max(theta, 1.0e-8_dp), 1.0_dp - 1.0e-8_dp)
        p(1) = log(max(omega, 1.0e-12_dp))
        p(2) = log((a / amax) / (1.0_dp - a / amax))
        p(3) = log(((g + 1.0_dp) / 3.0_dp) / (1.0_dp - (g + 1.0_dp) / 3.0_dp))
        p(4) = log(th / (1.0_dp - th))
    end subroutine midas_hyperbolic_asym_inv_transform

    subroutine midas_hyperbolic_weights(theta, weights, dweights)
        real(dp), intent(in) :: theta
        real(dp), intent(out) :: weights(midas_hyperbolic_m), dweights(midas_hyperbolic_m)
        real(dp) :: raw(midas_hyperbolic_m), draw(midas_hyperbolic_m)
        real(dp) :: ratio, sum_raw, sum_draw
        integer :: i

        raw(1) = theta
        draw(1) = 1.0_dp
        do i = 2, midas_hyperbolic_m
            ratio = (real(i - 1, dp) + theta) / real(i, dp)
            raw(i) = raw(i - 1) * ratio
            draw(i) = draw(i - 1) * ratio + raw(i - 1) / real(i, dp)
        end do
        sum_raw = sum(raw)
        sum_draw = sum(draw)
        weights = raw / sum_raw
        dweights = (draw*sum_raw - raw*sum_draw) / sum_raw**2
    end subroutine midas_hyperbolic_weights

    real(dp) function midas_hyperbolic_backcast(y) result(backcast)
        real(dp), intent(in) :: y(:)
        real(dp) :: weight, sum_weight
        integer :: i, tau

        tau = min(75, size(y))
        backcast = 0.0_dp
        sum_weight = 0.0_dp
        do i = 1, tau
            weight = 0.94_dp**real(i - 1, dp)
            backcast = backcast + weight*y(i)**2
            sum_weight = sum_weight + weight
        end do
        backcast = max(backcast / sum_weight, 1.0e-12_dp)
    end function midas_hyperbolic_backcast

    subroutine midas_hyperbolic_lag_terms(t, weights, dweights, backcast, x, dx_dtheta)
        integer, intent(in) :: t
        real(dp), intent(in) :: weights(midas_hyperbolic_m), dweights(midas_hyperbolic_m), backcast
        real(dp), intent(out) :: x, dx_dtheta
        real(dp) :: lag_sq
        integer :: i, idx

        x = 0.0_dp
        dx_dtheta = 0.0_dp
        do i = 1, midas_hyperbolic_m
            idx = t - i
            if (idx >= 1) then
                lag_sq = midas_hyperbolic_obs(idx)**2
            else
                lag_sq = backcast
            end if
            x = x + weights(i)*lag_sq
            dx_dtheta = dx_dtheta + dweights(i)*lag_sq
        end do
    end subroutine midas_hyperbolic_lag_terms

    subroutine midas_hyperbolic_asym_lag_terms(t, weights, dweights, backcast, x, xneg, dx_dtheta, dxneg_dtheta)
        integer, intent(in) :: t
        real(dp), intent(in) :: weights(midas_hyperbolic_m), dweights(midas_hyperbolic_m), backcast
        real(dp), intent(out) :: x, xneg, dx_dtheta, dxneg_dtheta
        real(dp) :: lag_sq, ind
        integer :: i, idx

        x = 0.0_dp
        xneg = 0.0_dp
        dx_dtheta = 0.0_dp
        dxneg_dtheta = 0.0_dp
        do i = 1, midas_hyperbolic_m
            idx = t - i
            if (idx >= 1) then
                lag_sq = midas_hyperbolic_obs(idx)**2
                ind = merge(1.0_dp, 0.0_dp, midas_hyperbolic_obs(idx) < 0.0_dp)
            else
                lag_sq = backcast
                ind = 0.5_dp
            end if
            x = x + weights(i)*lag_sq
            xneg = xneg + weights(i)*ind*lag_sq
            dx_dtheta = dx_dtheta + dweights(i)*lag_sq
            dxneg_dtheta = dxneg_dtheta + dweights(i)*ind*lag_sq
        end do
    end subroutine midas_hyperbolic_asym_lag_terms

    subroutine midas_hyperbolic_obj(p, np, f, g)
        integer, intent(in) :: np
        real(dp), intent(in) :: p(np)
        real(dp), intent(out) :: f, g(np)
        real(dp) :: omega, alpha, theta, weights(midas_hyperbolic_m), dweights(midas_hyperbolic_m)
        real(dp) :: backcast, x, dx_dtheta, h, y, factor
        real(dp) :: grad_om, grad_al, grad_th
        integer :: t

        call midas_hyperbolic_transform(p, omega, alpha, theta)
        call midas_hyperbolic_weights(theta, weights, dweights)
        backcast = midas_hyperbolic_backcast(midas_hyperbolic_obs)

        f = real(midas_hyperbolic_nobs, dp) * log_sqrt_2pi
        grad_om = 0.0_dp
        grad_al = 0.0_dp
        grad_th = 0.0_dp
        do t = 1, midas_hyperbolic_nobs
            call midas_hyperbolic_lag_terms(t, weights, dweights, backcast, x, dx_dtheta)
            h = max(omega + alpha*x, 1.0e-12_dp)
            y = midas_hyperbolic_obs(t)
            factor = 1.0_dp/h - y**2 / h**2
            f = f + 0.5_dp * (log(h) + y**2 / h)
            grad_om = grad_om + 0.5_dp * factor
            grad_al = grad_al + 0.5_dp * factor * x
            grad_th = grad_th + 0.5_dp * factor * alpha * dx_dtheta
        end do

        g(1) = grad_om * omega
        g(2) = grad_al * alpha * (1.0_dp - alpha)
        g(3) = grad_th * theta * (1.0_dp - theta)
        f = f / real(midas_hyperbolic_nobs, dp)
        g = g / real(midas_hyperbolic_nobs, dp)
    end subroutine midas_hyperbolic_obj

    subroutine midas_hyperbolic_asym_obj(p, np, f, g)
        integer, intent(in) :: np
        real(dp), intent(in) :: p(np)
        real(dp), intent(out) :: f, g(np)
        real(dp) :: omega, alpha, gamma, theta, weights(midas_hyperbolic_m), dweights(midas_hyperbolic_m)
        real(dp) :: backcast, x, xneg, dx_dtheta, dxneg_dtheta, h, y, factor
        real(dp) :: grad_om, grad_al, grad_ga, grad_th
        real(dp) :: amax, logistic_alpha, logistic_gamma, dalpha_dgamma, dgamma_dp
        integer :: t

        call midas_hyperbolic_asym_transform(p, omega, alpha, gamma, theta)
        call midas_hyperbolic_weights(theta, weights, dweights)
        backcast = midas_hyperbolic_backcast(midas_hyperbolic_obs)

        f = real(midas_hyperbolic_nobs, dp) * log_sqrt_2pi
        grad_om = 0.0_dp
        grad_al = 0.0_dp
        grad_ga = 0.0_dp
        grad_th = 0.0_dp
        do t = 1, midas_hyperbolic_nobs
            call midas_hyperbolic_asym_lag_terms(t, weights, dweights, backcast, x, xneg, dx_dtheta, dxneg_dtheta)
            h = max(omega + alpha*x + gamma*xneg, 1.0e-12_dp)
            y = midas_hyperbolic_obs(t)
            factor = 1.0_dp/h - y**2 / h**2
            f = f + 0.5_dp * (log(h) + y**2 / h)
            grad_om = grad_om + 0.5_dp * factor
            grad_al = grad_al + 0.5_dp * factor * x
            grad_ga = grad_ga + 0.5_dp * factor * xneg
            grad_th = grad_th + 0.5_dp * factor * (alpha*dx_dtheta + gamma*dxneg_dtheta)
        end do

        amax = min(1.0_dp - 0.5_dp*gamma, 1.0_dp)
        amax = max(amax, 1.0e-8_dp)
        logistic_alpha = alpha / amax
        logistic_gamma = (gamma + 1.0_dp) / 3.0_dp
        dgamma_dp = 3.0_dp * logistic_gamma * (1.0_dp - logistic_gamma)
        dalpha_dgamma = merge(-0.5_dp*logistic_alpha, 0.0_dp, gamma > 0.0_dp)

        g(1) = grad_om * omega
        g(2) = grad_al * amax * logistic_alpha * (1.0_dp - logistic_alpha)
        g(3) = (grad_ga + grad_al*dalpha_dgamma) * dgamma_dp
        g(4) = grad_th * theta * (1.0_dp - theta)
        f = f / real(midas_hyperbolic_nobs, dp)
        g = g / real(midas_hyperbolic_nobs, dp)
    end subroutine midas_hyperbolic_asym_obj

    subroutine aparch_transform(p, omega, alpha, gamma, beta, delta)
        real(dp), intent(in)  :: p(aparch_np)
        real(dp), intent(out) :: omega, alpha, gamma, beta, delta
        real(dp) :: e2, e4, s, u
        omega = exp(p(1))
        e2 = exp(p(2))
        e4 = exp(p(4))
        s = 1.0_dp + e2 + e4
        alpha = e2 / s
        beta = e4 / s
        gamma = 0.999_dp * tanh(p(3))
        u = 1.0_dp / (1.0_dp + exp(-p(5)))
        delta = 0.25_dp + 3.75_dp * u
    end subroutine aparch_transform

    subroutine aparch_inv_transform(omega, alpha, gamma, beta, delta, p)
        real(dp), intent(in)  :: omega, alpha, gamma, beta, delta
        real(dp), intent(out) :: p(aparch_np)
        real(dp) :: slack, g, d
        slack = max(1.0_dp - alpha - beta, 1.0e-8_dp)
        g = max(min(gamma / 0.999_dp, 0.999999_dp), -0.999999_dp)
        d = max(min(delta, 3.999999_dp), 0.250001_dp)
        p(1) = log(max(omega, 1.0e-12_dp))
        p(2) = log(max(alpha, 1.0e-12_dp) / slack)
        p(3) = 0.5_dp * log((1.0_dp + g) / (1.0_dp - g))
        p(4) = log(max(beta, 1.0e-12_dp) / slack)
        p(5) = log((d - 0.25_dp) / (4.0_dp - d))
    end subroutine aparch_inv_transform

    real(dp) function aparch_nll_value(p) result(f)
        real(dp), intent(in) :: p(aparch_np)
        real(dp) :: omega, alpha, gamma, beta, delta, sdel, h, zscale, term
        integer :: t

        call aparch_transform(p, omega, alpha, gamma, beta, delta)
        sdel = max(sum(abs(aparch_obs)**delta) / real(aparch_nobs, dp), 1.0e-12_dp)
        f = real(aparch_nobs, dp) * log_sqrt_2pi
        do t = 1, aparch_nobs
            h = max(sdel**(2.0_dp / delta), 1.0e-12_dp)
            f = f + 0.5_dp * (log(h) + aparch_obs(t)**2 / h)
            zscale = abs(aparch_obs(t)) - gamma*aparch_obs(t)
            term = max(zscale, 1.0e-12_dp)**delta
            sdel = omega + alpha*term + beta*sdel
            sdel = max(sdel, 1.0e-12_dp)
        end do
        f = f / real(aparch_nobs, dp)
    end function aparch_nll_value

    subroutine aparch_obj(p, np, f, g)
        integer,  intent(in)  :: np
        real(dp), intent(in)  :: p(np)
        real(dp), intent(out) :: f
        real(dp), intent(out) :: g(np)
        logical, parameter :: use_finite_difference_backup = .false.
        real(dp) :: pp(aparch_np), pm(aparch_np), fp, fm, step
        real(dp) :: omega, alpha, gamma, beta, delta
        real(dp) :: sdel, s_old, h, log_s, inv_s, y, ay, zscale, term, log_term
        real(dp) :: ds_dom, ds_dal, ds_dga, ds_dbe, ds_dde
        real(dp) :: ds_dom_old, ds_dal_old, ds_dga_old, ds_dbe_old, ds_dde_old
        real(dp) :: dh_dom, dh_dal, dh_dga, dh_dbe, dh_dde
        real(dp) :: grad_om, grad_al, grad_ga, grad_be, grad_de, factor
        real(dp) :: u, tg
        integer :: i, t

        if (use_finite_difference_backup) then
            f = aparch_nll_value(p)
            do i = 1, np
                pp = p
                pm = p
                step = 1.0e-5_dp * max(1.0_dp, abs(p(i)))
                pp(i) = pp(i) + step
                pm(i) = pm(i) - step
                fp = aparch_nll_value(pp)
                fm = aparch_nll_value(pm)
                g(i) = (fp - fm) / (2.0_dp * step)
            end do
            return
        end if

        call aparch_transform(p, omega, alpha, gamma, beta, delta)

        sdel = 0.0_dp
        ds_dde = 0.0_dp
        do t = 1, aparch_nobs
            ay = max(abs(aparch_obs(t)), 1.0e-12_dp)
            term = ay**delta
            sdel = sdel + term
            ds_dde = ds_dde + term*log(ay)
        end do
        sdel = max(sdel / real(aparch_nobs, dp), 1.0e-12_dp)
        ds_dom = 0.0_dp
        ds_dal = 0.0_dp
        ds_dga = 0.0_dp
        ds_dbe = 0.0_dp
        ds_dde = ds_dde / real(aparch_nobs, dp)

        f = real(aparch_nobs, dp) * log_sqrt_2pi
        grad_om = 0.0_dp
        grad_al = 0.0_dp
        grad_ga = 0.0_dp
        grad_be = 0.0_dp
        grad_de = 0.0_dp

        do t = 1, aparch_nobs
            sdel = max(sdel, 1.0e-12_dp)
            log_s = log(sdel)
            inv_s = 1.0_dp / sdel
            h = exp((2.0_dp / delta) * log_s)
            h = max(h, 1.0e-12_dp)

            dh_dom = h * (2.0_dp / delta) * ds_dom * inv_s
            dh_dal = h * (2.0_dp / delta) * ds_dal * inv_s
            dh_dga = h * (2.0_dp / delta) * ds_dga * inv_s
            dh_dbe = h * (2.0_dp / delta) * ds_dbe * inv_s
            dh_dde = h * ((2.0_dp / delta) * ds_dde * inv_s - (2.0_dp / delta**2) * log_s)

            y = aparch_obs(t)
            factor = 1.0_dp/h - y**2 / h**2
            f = f + 0.5_dp * (log(h) + y**2 / h)
            grad_om = grad_om + 0.5_dp * factor * dh_dom
            grad_al = grad_al + 0.5_dp * factor * dh_dal
            grad_ga = grad_ga + 0.5_dp * factor * dh_dga
            grad_be = grad_be + 0.5_dp * factor * dh_dbe
            grad_de = grad_de + 0.5_dp * factor * dh_dde

            s_old = sdel
            ds_dom_old = ds_dom
            ds_dal_old = ds_dal
            ds_dga_old = ds_dga
            ds_dbe_old = ds_dbe
            ds_dde_old = ds_dde

            zscale = max(abs(y) - gamma*y, 1.0e-12_dp)
            term = zscale**delta
            log_term = log(zscale)
            sdel = omega + alpha*term + beta*s_old
            ds_dom = 1.0_dp + beta*ds_dom_old
            ds_dal = term + beta*ds_dal_old
            ds_dga = -alpha*delta*y*zscale**(delta - 1.0_dp) + beta*ds_dga_old
            ds_dbe = s_old + beta*ds_dbe_old
            ds_dde = alpha*term*log_term + beta*ds_dde_old
        end do

        tg = tanh(p(3))
        u = 1.0_dp / (1.0_dp + exp(-p(5)))

        g(1) = grad_om * omega
        g(2) = alpha * (grad_al*(1.0_dp - alpha) - grad_be*beta)
        g(3) = grad_ga * 0.999_dp * (1.0_dp - tg**2)
        g(4) = beta * (-grad_al*alpha + grad_be*(1.0_dp - beta))
        g(5) = grad_de * 3.75_dp * u * (1.0_dp - u)

        f = f / real(aparch_nobs, dp)
        g = g / real(aparch_nobs, dp)
    end subroutine aparch_obj

    real(dp) function aparch_mean_variance(y, params)
        real(dp), intent(in) :: y(:)
        type(garch_params_t), intent(in) :: params
        real(dp) :: sdel, h, term
        integer :: t, n

        n = size(y)
        sdel = max(sum(abs(y)**params%theta) / real(n, dp), 1.0e-12_dp)
        aparch_mean_variance = 0.0_dp
        do t = 1, n
            h = max(sdel**(2.0_dp / params%theta), 1.0e-12_dp)
            aparch_mean_variance = aparch_mean_variance + h
            term = max(abs(y(t)) - params%gamma*y(t), 1.0e-12_dp)**params%theta
            sdel = params%omega + params%alpha*term + params%beta*sdel
            sdel = max(sdel, 1.0e-12_dp)
        end do
        aparch_mean_variance = aparch_mean_variance / real(n, dp)
    end function aparch_mean_variance

    subroutine garch_skew_kurt(y, params, skew, ekurt)
        real(dp), intent(in)  :: y(:)
        type(garch_params_t), intent(in) :: params
        real(dp), intent(out) :: skew, ekurt
        real(dp) :: h, z
        real(dp), allocatable :: zz(:)
        integer :: t, n

        n = size(y)
        h = params%omega / max(1.0_dp - symm_garch_persist(params), 1.0e-8_dp)
        allocate(zz(n))
        do t = 1, n
            h = max(h, 1.0e-12_dp)
            z = y(t) / sqrt(h)
            zz(t) = z
            h = params%omega + params%alpha*y(t)**2 + params%beta*h
        end do
        call moments(zz, n, skew, ekurt)
        deallocate(zz)
    end subroutine garch_skew_kurt

    subroutine nagarch_skew_kurt(y, params, skew, ekurt)
        real(dp), intent(in)  :: y(:)
        type(garch_params_t), intent(in) :: params
        real(dp), intent(out) :: skew, ekurt
        real(dp) :: h, sqrth, r, z
        real(dp), allocatable :: zz(:)
        integer :: t, n

        n = size(y)
        h = params%omega / max(1.0_dp - nagarch_persist(params), 1.0e-8_dp)
        allocate(zz(n))
        do t = 1, n
            h = max(h, 1.0e-12_dp)
            sqrth = sqrt(h)
            z = y(t) / sqrth
            zz(t) = z
            r = y(t) - params%theta*sqrth
            h = params%omega + params%alpha*r**2 + params%beta*h
        end do
        call moments(zz, n, skew, ekurt)
        deallocate(zz)
    end subroutine nagarch_skew_kurt

    subroutine qgarch_skew_kurt(y, params, skew, ekurt)
        real(dp), intent(in)  :: y(:)
        type(garch_params_t), intent(in) :: params
        real(dp), intent(out) :: skew, ekurt
        real(dp) :: h, z
        real(dp), allocatable :: zz(:)
        integer :: t, n

        n = size(y)
        h = qgarch_mean_variance(params)
        allocate(zz(n))
        do t = 1, n
            h = max(h, 1.0e-12_dp)
            z = y(t) / sqrt(h)
            zz(t) = z
            h = params%omega + params%alpha*(y(t) - params%theta)**2 + params%beta*h
        end do
        call moments(zz, n, skew, ekurt)
        deallocate(zz)
    end subroutine qgarch_skew_kurt

    subroutine rgarch_skew_kurt(y, range_var, params, skew, ekurt)
        real(dp), intent(in)  :: y(:), range_var(:)
        type(garch_params_t), intent(in) :: params
        real(dp), intent(out) :: skew, ekurt
        real(dp) :: h, z
        real(dp), allocatable :: zz(:)
        integer :: t, n

        n = size(y)
        h = rgarch_initial_variance(y, range_var, params)
        allocate(zz(n))
        do t = 1, n
            h = max(h, 1.0e-12_dp)
            z = y(t) / sqrt(h)
            zz(t) = z
            h = params%omega + params%alpha*range_var(t) + params%beta*h
        end do
        call moments(zz, n, skew, ekurt)
        deallocate(zz)
    end subroutine rgarch_skew_kurt

    subroutine carr_park_skew_kurt(y, range_var, params, skew, ekurt)
        real(dp), intent(in)  :: y(:), range_var(:)
        type(garch_params_t), intent(in) :: params
        real(dp), intent(out) :: skew, ekurt
        real(dp) :: h, z
        real(dp), allocatable :: zz(:)
        integer :: t, n

        n = size(y)
        if (size(range_var) /= n) then
            print '(A)', "garch_fit_mod: range_var length mismatch in carr_park_skew_kurt"
            error stop
        end if
        h = carr_park_initial_variance(range_var, params)
        allocate(zz(n))
        do t = 1, n
            h = max(h, 1.0e-12_dp)
            z = y(t) / sqrt(h)
            zz(t) = z
            h = params%omega + params%alpha*range_var(t) + params%beta*h
        end do
        call moments(zz, n, skew, ekurt)
        deallocate(zz)
    end subroutine carr_park_skew_kurt

    subroutine regarch1_skew_kurt(y, log_range, params, skew, ekurt)
        real(dp), intent(in)  :: y(:), log_range(:)
        type(garch_params_t), intent(in) :: params
        real(dp), intent(out) :: skew, ekurt
        real(dp) :: lv, sig, z, x_range
        real(dp), allocatable :: zz(:)
        integer :: t, n

        n = size(y)
        lv = regarch1_initial_log_volatility(log_range)
        allocate(zz(n))
        do t = 1, n
            sig = exp(lv)
            z = y(t) / max(sig, 1.0e-8_dp)
            zz(t) = z
            x_range = (log_range(t) - regarch_log_range_mean - lv) / regarch_log_range_sd
            lv = params%omega + params%beta*lv + params%alpha*x_range + params%gamma*z
        end do
        call moments(zz, n, skew, ekurt)
        deallocate(zz)
    end subroutine regarch1_skew_kurt

    subroutine regarch2_skew_kurt(y, log_range, params, skew, ekurt)
        real(dp), intent(in)  :: y(:), log_range(:)
        type(garch_params_t), intent(in) :: params
        real(dp), intent(out) :: skew, ekurt
        real(dp) :: lh, lq, sig, z, x_range, lh_next, lq_next
        real(dp), allocatable :: zz(:)
        integer :: t, n

        n = size(y)
        lh = regarch1_initial_log_volatility(log_range)
        lq = lh
        allocate(zz(n))
        do t = 1, n
            sig = exp(lh)
            z = y(t) / max(sig, 1.0e-8_dp)
            zz(t) = z
            x_range = (log_range(t) - regarch_log_range_mean - lh) / regarch_log_range_sd
            lh_next = params%beta*lh + (1.0_dp - params%beta)*lq + params%alpha*x_range + params%gamma*z
            lq_next = (1.0_dp - params%scale)*params%omega + params%scale*lq + &
                      params%theta*x_range + params%twist*z
            lh = lh_next
            lq = lq_next
        end do
        call moments(zz, n, skew, ekurt)
        deallocate(zz)
    end subroutine regarch2_skew_kurt

    subroutine rgarch_meas_skew_kurt(y, log_range, params, skew, ekurt)
        real(dp), intent(in)  :: y(:), log_range(:)
        type(garch_params_t), intent(in) :: params
        real(dp), intent(out) :: skew, ekurt
        real(dp) :: lv, sig, z, meas_hat, u
        real(dp), allocatable :: zz(:)
        integer :: t, n

        n = size(y)
        lv = log(max(sum(y**2) / real(n, dp), 1.0e-12_dp)) / 2.0_dp
        allocate(zz(n))
        do t = 1, n
            sig = exp(lv)
            z = y(t) / max(sig, 1.0e-8_dp)
            zz(t) = z
            meas_hat = params%theta + params%twist*lv + params%extra1*z + params%extra2*(z**2 - 1.0_dp)
            u = (log_range(t) - meas_hat) / max(params%scale, 1.0e-8_dp)
            lv = params%omega + params%beta*lv + params%alpha*u + params%gamma*z
        end do
        call moments(zz, n, skew, ekurt)
        deallocate(zz)
    end subroutine rgarch_meas_skew_kurt

    subroutine nagarch_range_skew_kurt(y, range_var, params, skew, ekurt)
        real(dp), intent(in)  :: y(:), range_var(:)
        type(garch_params_t), intent(in) :: params
        real(dp), intent(out) :: skew, ekurt
        real(dp) :: h, sqrth, r, z
        real(dp), allocatable :: zz(:)
        integer :: t, n

        n = size(y)
        h = nagarch_range_initial_variance(y, range_var, params)
        allocate(zz(n))
        do t = 1, n
            h = max(h, 1.0e-12_dp)
            sqrth = sqrt(h)
            z = y(t) / sqrth
            zz(t) = z
            r = y(t) - params%theta*sqrth
            h = params%omega + params%alpha*r**2 + params%gamma*range_var(t) + params%beta*h
        end do
        call moments(zz, n, skew, ekurt)
        deallocate(zz)
    end subroutine nagarch_range_skew_kurt

    subroutine fgarch_twist_range_skew_kurt(y, range_var, params, skew, ekurt)
        real(dp), intent(in)  :: y(:), range_var(:)
        type(garch_params_t), intent(in) :: params
        real(dp), intent(out) :: skew, ekurt
        real(dp) :: h, sqrth, z, q
        real(dp), allocatable :: zz(:)
        integer :: t, n

        n = size(y)
        h = fgarch_twist_range_initial_variance(y, range_var, params)
        allocate(zz(n))
        do t = 1, n
            h = max(h, 1.0e-12_dp)
            sqrth = sqrt(h)
            z = y(t) / sqrth
            zz(t) = z
            q = abs(z - params%theta) - params%twist*(z - params%theta)
            h = params%omega + params%alpha*h*q**2 + params%gamma*range_var(t) + params%beta*h
        end do
        call moments(zz, n, skew, ekurt)
        deallocate(zz)
    end subroutine fgarch_twist_range_skew_kurt

    subroutine gjr_skew_kurt(y, params, skew, ekurt)
        real(dp), intent(in)  :: y(:)
        type(garch_params_t), intent(in) :: params
        real(dp), intent(out) :: skew, ekurt
        real(dp) :: h, ind, z
        real(dp), allocatable :: zz(:)
        integer :: t, n

        n = size(y)
        h = params%omega / max(1.0_dp - gjr_persist(params), 1.0e-8_dp)
        allocate(zz(n))
        do t = 1, n
            h = max(h, 1.0e-12_dp)
            z = y(t) / sqrt(h)
            zz(t) = z
            ind = merge(1.0_dp, 0.0_dp, y(t) < 0.0_dp)
            h = params%omega + (params%alpha + params%gamma*ind)*y(t)**2 + params%beta*h
        end do
        call moments(zz, n, skew, ekurt)
        deallocate(zz)
    end subroutine gjr_skew_kurt

    subroutine egarch_skew_kurt(y, params, skew, ekurt)
        real(dp), intent(in)  :: y(:)
        type(garch_params_t), intent(in) :: params
        real(dp), intent(out) :: skew, ekurt
        real(dp) :: lh, h, z, c_eg
        real(dp), allocatable :: zz(:)
        integer :: t, n

        n = size(y)
        c_eg = sqrt(2.0_dp / 3.14159265358979323846_dp)
        lh = params%omega / (1.0_dp - params%beta)
        allocate(zz(n))
        do t = 1, n
            h = exp(lh)
            z = y(t) / sqrt(h)
            zz(t) = z
            lh = params%omega + params%beta*lh + params%alpha*(abs(z) - c_eg) + params%gamma*z
        end do
        call moments(zz, n, skew, ekurt)
        deallocate(zz)
    end subroutine egarch_skew_kurt

    subroutine aparch_skew_kurt(y, params, skew, ekurt)
        real(dp), intent(in)  :: y(:)
        type(garch_params_t), intent(in) :: params
        real(dp), intent(out) :: skew, ekurt
        real(dp) :: sdel, h, z, term
        real(dp), allocatable :: zz(:)
        integer :: t, n

        n = size(y)
        sdel = max(sum(abs(y)**params%theta) / real(n, dp), 1.0e-12_dp)
        allocate(zz(n))
        do t = 1, n
            h = max(sdel**(2.0_dp / params%theta), 1.0e-12_dp)
            z = y(t) / sqrt(h)
            zz(t) = z
            term = max(abs(y(t)) - params%gamma*y(t), 1.0e-12_dp)**params%theta
            sdel = params%omega + params%alpha*term + params%beta*sdel
            sdel = max(sdel, 1.0e-12_dp)
        end do
        call moments(zz, n, skew, ekurt)
        deallocate(zz)
    end subroutine aparch_skew_kurt

    subroutine harch_skew_kurt(y, params, skew, ekurt)
        real(dp), intent(in)  :: y(:)
        type(garch_params_t), intent(in) :: params
        real(dp), intent(out) :: skew, ekurt
        real(dp), allocatable :: zz(:)
        real(dp) :: backcast, h, z
        integer :: t, n

        n = size(y)
        call harch_set_data(y)
        backcast = max(sum(y**2) / real(n, dp), 1.0e-12_dp)
        allocate(zz(n))
        do t = 1, n
            h = params%omega + params%alpha*harch_block(t, 1, backcast) + &
                params%gamma*harch_block(t, 5, backcast) + params%beta*harch_block(t, 22, backcast)
            h = max(h, 1.0e-12_dp)
            z = y(t) / sqrt(h)
            zz(t) = z
        end do
        call moments(zz, n, skew, ekurt)
        deallocate(zz)
    end subroutine harch_skew_kurt

    subroutine figarch_skew_kurt(y, params, skew, ekurt)
        real(dp), intent(in)  :: y(:)
        type(garch_params_t), intent(in) :: params
        real(dp), intent(out) :: skew, ekurt
        real(dp), allocatable :: variance(:), zz(:)
        integer :: t, n

        n = size(y)
        allocate(variance(n), zz(n))
        call figarch_variance(y, params, variance)
        do t = 1, n
            zz(t) = y(t) / sqrt(max(variance(t), 1.0e-12_dp))
        end do
        call moments(zz, n, skew, ekurt)
        deallocate(variance, zz)
    end subroutine figarch_skew_kurt

    subroutine fi_nagarch_skew_kurt(y, params, skew, ekurt)
        real(dp), intent(in)  :: y(:)
        type(garch_params_t), intent(in) :: params
        real(dp), intent(out) :: skew, ekurt
        real(dp), allocatable :: variance(:), zz(:)
        integer :: t, n

        n = size(y)
        allocate(variance(n), zz(n))
        call fi_nagarch_variance(y, params, variance)
        do t = 1, n
            zz(t) = y(t) / sqrt(max(variance(t), 1.0e-12_dp))
        end do
        call moments(zz, n, skew, ekurt)
        deallocate(variance, zz)
    end subroutine fi_nagarch_skew_kurt

    subroutine riskmetrics2006_skew_kurt(y, skew, ekurt)
        real(dp), intent(in)  :: y(:)
        real(dp), intent(out) :: skew, ekurt
        real(dp), allocatable :: variance(:), zz(:)
        integer :: t, n

        n = size(y)
        allocate(variance(n), zz(n))
        call riskmetrics2006_variance(y, variance)
        do t = 1, n
            zz(t) = y(t) / sqrt(max(variance(t), 1.0e-12_dp))
        end do
        call moments(zz, n, skew, ekurt)
        deallocate(variance, zz)
    end subroutine riskmetrics2006_skew_kurt

    subroutine midas_hyperbolic_skew_kurt(y, params, skew, ekurt)
        real(dp), intent(in)  :: y(:)
        type(garch_params_t), intent(in) :: params
        real(dp), intent(out) :: skew, ekurt
        real(dp) :: weights(midas_hyperbolic_m), dweights(midas_hyperbolic_m)
        real(dp) :: backcast, x, dx_dtheta, h
        real(dp), allocatable :: zz(:)
        integer :: t, n

        n = size(y)
        call midas_hyperbolic_set_data(y)
        call midas_hyperbolic_weights(params%theta, weights, dweights)
        backcast = midas_hyperbolic_backcast(y)
        allocate(zz(n))
        do t = 1, n
            call midas_hyperbolic_lag_terms(t, weights, dweights, backcast, x, dx_dtheta)
            h = max(params%omega + params%alpha*x, 1.0e-12_dp)
            zz(t) = y(t) / sqrt(h)
        end do
        call moments(zz, n, skew, ekurt)
        deallocate(zz)
    end subroutine midas_hyperbolic_skew_kurt

    subroutine midas_hyperbolic_asym_skew_kurt(y, params, skew, ekurt)
        real(dp), intent(in)  :: y(:)
        type(garch_params_t), intent(in) :: params
        real(dp), intent(out) :: skew, ekurt
        real(dp) :: weights(midas_hyperbolic_m), dweights(midas_hyperbolic_m)
        real(dp) :: backcast, h, lag_sq, ind
        real(dp), allocatable :: zz(:)
        integer :: t, i, idx, n

        n = size(y)
        call midas_hyperbolic_set_data(y)
        call midas_hyperbolic_weights(params%theta, weights, dweights)
        backcast = midas_hyperbolic_backcast(y)
        allocate(zz(n))
        do t = 1, n
            h = params%omega
            do i = 1, midas_hyperbolic_m
                idx = t - i
                if (idx >= 1) then
                    lag_sq = y(idx)**2
                    ind = merge(1.0_dp, 0.0_dp, y(idx) < 0.0_dp)
                else
                    lag_sq = backcast
                    ind = 0.5_dp
                end if
                h = h + (params%alpha + params%gamma*ind)*weights(i)*lag_sq
            end do
            h = max(h, 1.0e-12_dp)
            zz(t) = y(t) / sqrt(h)
        end do
        call moments(zz, n, skew, ekurt)
        deallocate(zz)
    end subroutine midas_hyperbolic_asym_skew_kurt

    subroutine ewma_skew_kurt(y, params, skew, ekurt)
        real(dp), intent(in)  :: y(:)
        type(garch_params_t), intent(in) :: params
        real(dp), intent(out) :: skew, ekurt
        real(dp) :: h, z
        real(dp), allocatable :: zz(:)
        integer :: t, n

        n = size(y)
        h = sample_variance(y)
        allocate(zz(n))
        do t = 1, n
            h = max(h, 1.0e-12_dp)
            z = y(t) / sqrt(h)
            zz(t) = z
            h = params%beta*h + params%alpha*y(t)**2
        end do
        call moments(zz, n, skew, ekurt)
        deallocate(zz)
    end subroutine ewma_skew_kurt

    subroutine aewma_nag_skew_kurt(y, params, skew, ekurt)
        real(dp), intent(in)  :: y(:)
        type(garch_params_t), intent(in) :: params
        real(dp), intent(out) :: skew, ekurt
        real(dp) :: h, sqrth, z, nic
        real(dp), allocatable :: zz(:)
        integer :: t, n

        n = size(y)
        h = sample_variance(y)
        allocate(zz(n))
        do t = 1, n
            h = max(h, 1.0e-12_dp)
            sqrth = sqrt(h)
            z = y(t) / sqrth
            zz(t) = z
            nic = (z - params%theta)**2
            h = params%beta*h + params%alpha*params%scale*h*nic
        end do
        call moments(zz, n, skew, ekurt)
        deallocate(zz)
    end subroutine aewma_nag_skew_kurt

    subroutine aewma_twist_skew_kurt(y, params, skew, ekurt)
        real(dp), intent(in)  :: y(:)
        type(garch_params_t), intent(in) :: params
        real(dp), intent(out) :: skew, ekurt
        real(dp) :: h, sqrth, z, q, nic
        real(dp), allocatable :: zz(:)
        integer :: t, n

        n = size(y)
        h = sample_variance(y)
        allocate(zz(n))
        do t = 1, n
            h = max(h, 1.0e-12_dp)
            sqrth = sqrt(h)
            z = y(t) / sqrth
            zz(t) = z
            q = abs(z - params%theta) - params%twist*(z - params%theta)
            nic = q**2
            h = params%beta*h + params%alpha*params%scale*h*nic
        end do
        call moments(zz, n, skew, ekurt)
        deallocate(zz)
    end subroutine aewma_twist_skew_kurt

    real(dp) function fgarch_twist_moment(theta, twist)
        real(dp), intent(in) :: theta, twist
        real(dp) :: Phi, ph, A, B, one_minus, one_plus

        Phi = 0.5_dp * (1.0_dp + erf(theta / sqrt2))
        ph  = exp(-0.5_dp * theta**2 - log_sqrt_2pi)
        A = (1.0_dp + theta**2) * (1.0_dp - Phi) - theta * ph
        B = (1.0_dp + theta**2) * Phi + theta * ph
        one_minus = 1.0_dp - twist
        one_plus  = 1.0_dp + twist
        fgarch_twist_moment = one_minus**2 * A + one_plus**2 * B
    end function fgarch_twist_moment

    real(dp) function symm_garch_persist(params)
        type(garch_params_t), intent(in) :: params
        symm_garch_persist = params%alpha + params%beta
    end function symm_garch_persist

    real(dp) function qgarch_persist(params)
        type(garch_params_t), intent(in) :: params
        qgarch_persist = params%alpha + params%beta
    end function qgarch_persist

    real(dp) function qgarch_mean_variance(params)
        type(garch_params_t), intent(in) :: params
        qgarch_mean_variance = qgarch_initial_variance(params%omega, params%alpha, params%beta, params%theta)
    end function qgarch_mean_variance

    real(dp) function nagarch_persist(params)
        type(garch_params_t), intent(in) :: params
        nagarch_persist = params%alpha * (1.0_dp + params%theta**2) + params%beta
    end function nagarch_persist

    real(dp) function rgarch_persist(params)
        type(garch_params_t), intent(in) :: params
        rgarch_persist = params%alpha + params%beta
    end function rgarch_persist

    real(dp) function carr_park_persist(params)
        type(garch_params_t), intent(in) :: params
        carr_park_persist = params%alpha + params%beta
    end function carr_park_persist

    real(dp) function regarch1_persist(params)
        type(garch_params_t), intent(in) :: params
        regarch1_persist = params%beta
    end function regarch1_persist

    real(dp) function regarch2_persist(params)
        type(garch_params_t), intent(in) :: params
        regarch2_persist = max(params%beta, params%scale)
    end function regarch2_persist

    real(dp) function rgarch_meas_persist(params)
        type(garch_params_t), intent(in) :: params
        rgarch_meas_persist = params%beta
    end function rgarch_meas_persist

    real(dp) function gjr_persist(params)
        type(garch_params_t), intent(in) :: params
        gjr_persist = params%alpha + 0.5_dp*params%gamma + params%beta
    end function gjr_persist

    real(dp) function egarch_persist(params)
        type(garch_params_t), intent(in) :: params
        egarch_persist = params%beta
    end function egarch_persist

    real(dp) function aparch_persist(params)
        type(garch_params_t), intent(in) :: params
        aparch_persist = params%alpha + params%beta
    end function aparch_persist

    real(dp) function harch_persist(params)
        type(garch_params_t), intent(in) :: params
        harch_persist = params%alpha + params%gamma + params%beta
    end function harch_persist

    real(dp) function figarch_persist(params)
        type(garch_params_t), intent(in) :: params
        real(dp), allocatable :: lambda(:)

        allocate(lambda(figarch_trunc_lag))
        call figarch_weights(params%alpha, params%theta, params%beta, lambda)
        figarch_persist = sum(lambda)
        deallocate(lambda)
    end function figarch_persist

    real(dp) function fi_nagarch_persist(params)
        type(garch_params_t), intent(in) :: params
        fi_nagarch_persist = figarch_persist(params)
    end function fi_nagarch_persist

    real(dp) function riskmetrics2006_persist()
        riskmetrics2006_persist = 1.0_dp
    end function riskmetrics2006_persist

    real(dp) function midas_hyperbolic_persist(params)
        type(garch_params_t), intent(in) :: params
        midas_hyperbolic_persist = params%alpha
    end function midas_hyperbolic_persist

    real(dp) function midas_hyperbolic_asym_persist(params)
        type(garch_params_t), intent(in) :: params
        midas_hyperbolic_asym_persist = params%alpha + 0.5_dp*params%gamma
    end function midas_hyperbolic_asym_persist

    real(dp) function fgarch_twist_persist(params)
        type(garch_params_t), intent(in) :: params
        fgarch_twist_persist = params%alpha * fgarch_twist_moment(params%theta, params%twist) + params%beta
    end function fgarch_twist_persist

    real(dp) function ewma_persist(params)
        type(garch_params_t), intent(in) :: params
        ewma_persist = params%alpha + params%beta
    end function ewma_persist

    real(dp) function aewma_nag_persist(params)
        type(garch_params_t), intent(in) :: params
        aewma_nag_persist = params%beta + params%alpha*params%scale*(1.0_dp + params%theta**2)
    end function aewma_nag_persist

    real(dp) function aewma_twist_persist(params)
        type(garch_params_t), intent(in) :: params
        aewma_twist_persist = params%beta + params%alpha*params%scale* &
            fgarch_twist_moment(params%theta, params%twist)
    end function aewma_twist_persist

    subroutine ewma_set_data(y)
        real(dp), intent(in) :: y(:)

        ewma_nobs = size(y)
        if (allocated(ewma_obs)) deallocate(ewma_obs)
        allocate(ewma_obs(ewma_nobs))
        ewma_obs = y
    end subroutine ewma_set_data

    subroutine ewma_obj(p, np, f, g)
        integer,  intent(in)  :: np
        real(dp), intent(in)  :: p(np)
        real(dp), intent(out) :: f, g(np)
        real(dp) :: lambda, dlambda_dp, h, dh, yt2, fsum, gsum, dfdh
        integer :: t

        call ewma_transform(p, lambda)
        dlambda_dp = lambda * (1.0_dp - lambda)
        h = sample_variance(ewma_obs)
        dh = 0.0_dp
        fsum = 0.0_dp
        gsum = 0.0_dp

        do t = 1, ewma_nobs
            h = max(h, 1.0e-12_dp)
            yt2 = ewma_obs(t)**2
            fsum = fsum + log_sqrt_2pi + 0.5_dp*(log(h) + yt2/h)
            dfdh = 0.5_dp*(1.0_dp - yt2/h) / h
            gsum = gsum + dfdh*dh
            dh = h + lambda*dh - yt2
            h = lambda*h + (1.0_dp - lambda)*yt2
        end do

        f = fsum / real(ewma_nobs, dp)
        g(1) = gsum * dlambda_dp / real(ewma_nobs, dp)
    end subroutine ewma_obj

    subroutine ewma_transform(p, lambda)
        real(dp), intent(in)  :: p(ewma_np)
        real(dp), intent(out) :: lambda

        lambda = 1.0_dp / (1.0_dp + exp(-p(1)))
        lambda = min(max(lambda, 1.0e-6_dp), 1.0_dp - 1.0e-6_dp)
    end subroutine ewma_transform

    subroutine ewma_inv_transform(lambda, p)
        real(dp), intent(in)  :: lambda
        real(dp), intent(out) :: p(ewma_np)
        real(dp) :: lam

        lam = min(max(lambda, 1.0e-6_dp), 1.0_dp - 1.0e-6_dp)
        p(1) = log(lam / (1.0_dp - lam))
    end subroutine ewma_inv_transform

    subroutine rgarch_set_data(y, range_var)
        real(dp), intent(in) :: y(:), range_var(:)

        rgarch_nobs = size(y)
        if (size(range_var) /= rgarch_nobs) then
            print '(A)', "garch_fit_mod: range_var length mismatch in fit_rgarch"
            error stop
        end if
        if (allocated(rgarch_obs)) deallocate(rgarch_obs)
        if (allocated(rgarch_x)) deallocate(rgarch_x)
        allocate(rgarch_obs(rgarch_nobs), rgarch_x(rgarch_nobs))
        rgarch_obs = y
        rgarch_x = range_var
    end subroutine rgarch_set_data

    subroutine rgarch_obj(p, np, f, g)
        integer,  intent(in)  :: np
        real(dp), intent(in)  :: p(np)
        real(dp), intent(out) :: f, g(np)
        real(dp) :: omega, alpha, beta, h, factor
        real(dp) :: dh_dom, dh_dal, dh_dbe
        real(dp) :: grad_om, grad_al, grad_be
        integer :: t

        call rgarch_transform(p, omega, alpha, beta)
        h = rgarch_initial_variance_raw(rgarch_obs, rgarch_x, omega, alpha, beta)
        dh_dom = 1.0_dp / max(1.0_dp - beta, 1.0e-8_dp)
        dh_dal = sum(rgarch_x) / real(rgarch_nobs, dp) / max(1.0_dp - beta, 1.0e-8_dp)
        dh_dbe = (omega + alpha*sum(rgarch_x)/real(rgarch_nobs, dp)) / max(1.0_dp - beta, 1.0e-8_dp)**2

        f = real(rgarch_nobs, dp) * log_sqrt_2pi
        grad_om = 0.0_dp
        grad_al = 0.0_dp
        grad_be = 0.0_dp
        do t = 1, rgarch_nobs
            if (h <= 0.0_dp .or. h /= h .or. h > 1.0e100_dp) then
                f = huge(1.0_dp)
                g = 0.0_dp
                return
            end if
            h = max(h, 1.0e-12_dp)
            factor = 1.0_dp/h - rgarch_obs(t)**2/h**2
            f = f + 0.5_dp*(log(h) + rgarch_obs(t)**2/h)
            grad_om = grad_om + 0.5_dp*factor*dh_dom
            grad_al = grad_al + 0.5_dp*factor*dh_dal
            grad_be = grad_be + 0.5_dp*factor*dh_dbe
            dh_dom = 1.0_dp + beta*dh_dom
            dh_dal = rgarch_x(t) + beta*dh_dal
            dh_dbe = h + beta*dh_dbe
            h = omega + alpha*rgarch_x(t) + beta*h
        end do

        g(1) = grad_om * omega
        g(2) = grad_al*alpha*(1.0_dp - alpha) - grad_be*beta*alpha
        g(3) = -grad_al*alpha*beta + grad_be*beta*(1.0_dp - beta)
        f = f / real(rgarch_nobs, dp)
        g = g / real(rgarch_nobs, dp)
    end subroutine rgarch_obj

    subroutine rgarch_transform(p, omega, alpha, beta)
        real(dp), intent(in)  :: p(:)
        real(dp), intent(out) :: omega, alpha, beta
        real(dp) :: e2, e3, s

        omega = exp(p(1))
        e2 = exp(p(2))
        e3 = exp(p(3))
        s = 1.0_dp + e2 + e3
        alpha = e2 / s
        beta = e3 / s
    end subroutine rgarch_transform

    subroutine rgarch_inv_transform(omega, alpha, beta, p)
        real(dp), intent(in)  :: omega, alpha, beta
        real(dp), intent(out) :: p(rgarch_np)
        real(dp) :: slack

        slack = max(1.0_dp - alpha - beta, 1.0e-8_dp)
        p(1) = log(max(omega, 1.0e-12_dp))
        p(2) = log(max(alpha, 1.0e-12_dp) / slack)
        p(3) = log(max(beta, 1.0e-12_dp) / slack)
    end subroutine rgarch_inv_transform

    real(dp) function rgarch_initial_variance(y, range_var, params)
        real(dp), intent(in) :: y(:), range_var(:)
        type(garch_params_t), intent(in) :: params

        rgarch_initial_variance = rgarch_initial_variance_raw(y, range_var, &
            params%omega, params%alpha, params%beta)
    end function rgarch_initial_variance

    real(dp) function rgarch_initial_variance_raw(y, range_var, omega, alpha, beta)
        real(dp), intent(in) :: y(:), range_var(:), omega, alpha, beta

        rgarch_initial_variance_raw = max((omega + alpha*sum(range_var)/real(size(range_var), dp)) / &
            max(1.0_dp - beta, 1.0e-8_dp), sum(y**2)/real(size(y), dp), 1.0e-12_dp)
    end function rgarch_initial_variance_raw

    subroutine carr_park_set_data(range_var)
        real(dp), intent(in) :: range_var(:)

        carr_park_nobs = size(range_var)
        if (allocated(carr_park_x)) deallocate(carr_park_x)
        allocate(carr_park_x(carr_park_nobs))
        carr_park_x = max(range_var, 1.0e-12_dp)
    end subroutine carr_park_set_data

    subroutine carr_park_obj(p, np, f, g)
        integer,  intent(in)  :: np
        real(dp), intent(in)  :: p(np)
        real(dp), intent(out) :: f, g(np)
        real(dp) :: omega, alpha, beta, h, x, factor
        real(dp) :: dh_dom, dh_dal, dh_dbe
        real(dp) :: grad_om, grad_al, grad_be
        integer :: t

        call rgarch_transform(p, omega, alpha, beta)
        h = carr_park_initial_variance_raw(carr_park_x, omega, alpha, beta)
        dh_dom = 1.0_dp / max(1.0_dp - beta, 1.0e-8_dp)
        dh_dal = sum(carr_park_x) / real(carr_park_nobs, dp) / max(1.0_dp - beta, 1.0e-8_dp)
        dh_dbe = (omega + alpha*sum(carr_park_x)/real(carr_park_nobs, dp)) / &
                 max(1.0_dp - beta, 1.0e-8_dp)**2

        f = 0.0_dp
        grad_om = 0.0_dp
        grad_al = 0.0_dp
        grad_be = 0.0_dp
        do t = 1, carr_park_nobs
            if (h <= 0.0_dp .or. h /= h .or. h > 1.0e100_dp) then
                f = huge(1.0_dp)
                g = 0.0_dp
                return
            end if
            h = max(h, 1.0e-12_dp)
            x = carr_park_x(t)
            factor = 1.0_dp/h - x/h**2
            f = f + log(h) + x/h
            grad_om = grad_om + factor*dh_dom
            grad_al = grad_al + factor*dh_dal
            grad_be = grad_be + factor*dh_dbe
            dh_dom = 1.0_dp + beta*dh_dom
            dh_dal = x + beta*dh_dal
            dh_dbe = h + beta*dh_dbe
            h = omega + alpha*x + beta*h
        end do

        g(1) = grad_om * omega
        g(2) = grad_al*alpha*(1.0_dp - alpha) - grad_be*beta*alpha
        g(3) = -grad_al*alpha*beta + grad_be*beta*(1.0_dp - beta)
        f = f / real(carr_park_nobs, dp)
        g = g / real(carr_park_nobs, dp)
    end subroutine carr_park_obj

    real(dp) function carr_park_initial_variance(range_var, params)
        real(dp), intent(in) :: range_var(:)
        type(garch_params_t), intent(in) :: params

        carr_park_initial_variance = carr_park_initial_variance_raw(range_var, &
            params%omega, params%alpha, params%beta)
    end function carr_park_initial_variance

    real(dp) function carr_park_initial_variance_raw(range_var, omega, alpha, beta)
        real(dp), intent(in) :: range_var(:), omega, alpha, beta
        real(dp) :: xbar

        xbar = sum(max(range_var, 1.0e-12_dp)) / real(size(range_var), dp)
        carr_park_initial_variance_raw = max((omega + alpha*xbar) / &
            max(1.0_dp - beta, 1.0e-8_dp), xbar, 1.0e-12_dp)
    end function carr_park_initial_variance_raw

    subroutine regarch1_set_data(y, log_range)
        real(dp), intent(in) :: y(:), log_range(:)

        regarch_nobs = size(y)
        if (size(log_range) /= regarch_nobs) then
            print '(A)', "garch_fit_mod: log_range length mismatch in fit_regarch1"
            error stop
        end if
        if (allocated(regarch_obs)) deallocate(regarch_obs)
        if (allocated(regarch_log_range)) deallocate(regarch_log_range)
        allocate(regarch_obs(regarch_nobs), regarch_log_range(regarch_nobs))
        regarch_obs = y
        regarch_log_range = log_range
    end subroutine regarch1_set_data

    subroutine regarch1_obj(p, np, f, g)
        integer,  intent(in)  :: np
        real(dp), intent(in)  :: p(np)
        real(dp), intent(out) :: f, g(np)

        real(dp) :: omega, alpha, gamma, beta, lv, sig, z, x_range, resid, factor, kappa
        real(dp) :: dlv_dom, dlv_dal, dlv_dga, dlv_dbe
        real(dp) :: grad_om, grad_al, grad_ga, grad_be
        integer :: t

        call regarch1_transform(p, omega, alpha, gamma, beta)
        lv = regarch1_initial_log_volatility(regarch_log_range)
        dlv_dom = 0.0_dp
        dlv_dal = 0.0_dp
        dlv_dga = 0.0_dp
        dlv_dbe = 0.0_dp

        f = real(regarch_nobs, dp) * (log_sqrt_2pi + log(regarch_log_range_sd))
        grad_om = 0.0_dp
        grad_al = 0.0_dp
        grad_ga = 0.0_dp
        grad_be = 0.0_dp
        do t = 1, regarch_nobs
            if (lv /= lv .or. abs(lv) > 100.0_dp) then
                f = huge(1.0_dp)
                g = 0.0_dp
                return
            end if
            sig = exp(lv)
            z = regarch_obs(t) / max(sig, 1.0e-8_dp)
            x_range = (regarch_log_range(t) - regarch_log_range_mean - lv) / regarch_log_range_sd
            resid = regarch_log_range(t) - regarch_log_range_mean - lv
            factor = -resid / regarch_log_range_sd**2
            f = f + 0.5_dp * (resid / regarch_log_range_sd)**2
            grad_om = grad_om + factor * dlv_dom
            grad_al = grad_al + factor * dlv_dal
            grad_ga = grad_ga + factor * dlv_dga
            grad_be = grad_be + factor * dlv_dbe

            kappa = beta - alpha / regarch_log_range_sd - gamma*z
            dlv_dom = 1.0_dp + kappa*dlv_dom
            dlv_dal = x_range + kappa*dlv_dal
            dlv_dga = z + kappa*dlv_dga
            dlv_dbe = lv + kappa*dlv_dbe
            lv = omega + beta*lv + alpha*x_range + gamma*z
        end do

        g(1) = grad_om
        g(2) = grad_al
        g(3) = grad_ga
        g(4) = grad_be * beta * (1.0_dp - beta)
        f = f / real(regarch_nobs, dp)
        g = g / real(regarch_nobs, dp)
    end subroutine regarch1_obj

    subroutine regarch1_transform(p, omega, alpha, gamma, beta)
        real(dp), intent(in)  :: p(:)
        real(dp), intent(out) :: omega, alpha, gamma, beta

        omega = p(1)
        alpha = p(2)
        gamma = p(3)
        beta = 1.0_dp / (1.0_dp + exp(-p(4)))
        beta = min(max(beta, 1.0e-6_dp), 1.0_dp - 1.0e-6_dp)
    end subroutine regarch1_transform

    subroutine regarch1_inv_transform(omega, alpha, gamma, beta, p)
        real(dp), intent(in)  :: omega, alpha, gamma, beta
        real(dp), intent(out) :: p(regarch1_np)
        real(dp) :: b

        b = min(max(beta, 1.0e-6_dp), 1.0_dp - 1.0e-6_dp)
        p(1) = omega
        p(2) = alpha
        p(3) = gamma
        p(4) = log(b / (1.0_dp - b))
    end subroutine regarch1_inv_transform

    real(dp) function regarch1_initial_log_volatility(log_range)
        real(dp), intent(in) :: log_range(:)

        regarch1_initial_log_volatility = sum(log_range - regarch_log_range_mean) / real(size(log_range), dp)
    end function regarch1_initial_log_volatility

    subroutine regarch2_obj(p, np, f, g)
        integer,  intent(in)  :: np
        real(dp), intent(in)  :: p(np)
        real(dp), intent(out) :: f, g(np)
        real(dp) :: mu, alpha_h, gamma_h, beta_h, alpha_q, gamma_q, beta_q
        real(dp) :: lh, lq, sig, z, x_range, resid, factor, kappa_h, kappa_q
        real(dp) :: dlh(np), dlq(np), dlh_next(np), dlq_next(np), grad(np)
        integer :: t

        call regarch2_transform(p, mu, alpha_h, gamma_h, beta_h, alpha_q, gamma_q, beta_q)
        lh = regarch1_initial_log_volatility(regarch_log_range)
        lq = lh
        dlh = 0.0_dp
        dlq = 0.0_dp
        grad = 0.0_dp
        f = real(regarch_nobs, dp) * (log_sqrt_2pi + log(regarch_log_range_sd))
        do t = 1, regarch_nobs
            if (lh /= lh .or. lq /= lq .or. abs(lh) > 100.0_dp .or. abs(lq) > 100.0_dp) then
                f = huge(1.0_dp)
                g = 0.0_dp
                return
            end if
            sig = exp(lh)
            z = regarch_obs(t) / max(sig, 1.0e-8_dp)
            x_range = (regarch_log_range(t) - regarch_log_range_mean - lh) / regarch_log_range_sd
            resid = regarch_log_range(t) - regarch_log_range_mean - lh
            factor = -resid / regarch_log_range_sd**2
            f = f + 0.5_dp * (resid / regarch_log_range_sd)**2
            grad = grad + factor*dlh

            kappa_h = beta_h - alpha_h/regarch_log_range_sd - gamma_h*z
            kappa_q = -alpha_q/regarch_log_range_sd - gamma_q*z
            dlh_next = kappa_h*dlh + (1.0_dp - beta_h)*dlq
            dlq_next = kappa_q*dlh + beta_q*dlq
            dlh_next(2) = dlh_next(2) + x_range
            dlh_next(3) = dlh_next(3) + z
            dlh_next(4) = dlh_next(4) + lh - lq
            dlq_next(1) = dlq_next(1) + 1.0_dp - beta_q
            dlq_next(5) = dlq_next(5) + x_range
            dlq_next(6) = dlq_next(6) + z
            dlq_next(7) = dlq_next(7) + lq - mu

            lh = beta_h*lh + (1.0_dp - beta_h)*lq + alpha_h*x_range + gamma_h*z
            lq = (1.0_dp - beta_q)*mu + beta_q*lq + alpha_q*x_range + gamma_q*z
            dlh = dlh_next
            dlq = dlq_next
        end do

        g = grad
        g(4) = g(4) * beta_h * (1.0_dp - beta_h)
        g(7) = g(7) * beta_q * (1.0_dp - beta_q)
        f = f / real(regarch_nobs, dp)
        g = g / real(regarch_nobs, dp)
    end subroutine regarch2_obj

    subroutine regarch2_transform(p, mu, alpha_h, gamma_h, beta_h, alpha_q, gamma_q, beta_q)
        real(dp), intent(in)  :: p(:)
        real(dp), intent(out) :: mu, alpha_h, gamma_h, beta_h, alpha_q, gamma_q, beta_q

        mu = p(1)
        alpha_h = p(2)
        gamma_h = p(3)
        beta_h = 1.0_dp / (1.0_dp + exp(-p(4)))
        beta_h = min(max(beta_h, 1.0e-6_dp), 1.0_dp - 1.0e-6_dp)
        alpha_q = p(5)
        gamma_q = p(6)
        beta_q = 1.0_dp / (1.0_dp + exp(-p(7)))
        beta_q = min(max(beta_q, 1.0e-6_dp), 1.0_dp - 1.0e-6_dp)
    end subroutine regarch2_transform

    subroutine regarch2_inv_transform(mu, alpha_h, gamma_h, beta_h, alpha_q, gamma_q, beta_q, p)
        real(dp), intent(in)  :: mu, alpha_h, gamma_h, beta_h, alpha_q, gamma_q, beta_q
        real(dp), intent(out) :: p(regarch2_np)
        real(dp) :: bh, bq

        bh = min(max(beta_h, 1.0e-6_dp), 1.0_dp - 1.0e-6_dp)
        bq = min(max(beta_q, 1.0e-6_dp), 1.0_dp - 1.0e-6_dp)
        p(1) = mu
        p(2) = alpha_h
        p(3) = gamma_h
        p(4) = log(bh / (1.0_dp - bh))
        p(5) = alpha_q
        p(6) = gamma_q
        p(7) = log(bq / (1.0_dp - bq))
    end subroutine regarch2_inv_transform

    subroutine rgarch_meas_obj(p, np, f, g)
        integer,  intent(in)  :: np
        real(dp), intent(in)  :: p(np)
        real(dp), intent(out) :: f, g(np)
        real(dp) :: omega, alpha, gamma, beta, xi, phi, tau1, tau2, sigma_u
        real(dp) :: lv, sig, z, meas_hat, resid, u, z2, mcoeff, factor_lv
        real(dp) :: d(np), dz(np), dm(np), du(np), dnext(np), grad(np)
        integer :: t

        call rgarch_meas_transform(p, omega, alpha, gamma, beta, xi, phi, tau1, tau2, sigma_u)
        lv = log(max(sum(regarch_obs**2) / real(regarch_nobs, dp), 1.0e-12_dp)) / 2.0_dp
        d = 0.0_dp
        grad = 0.0_dp
        f = 0.0_dp

        do t = 1, regarch_nobs
            if (lv /= lv .or. abs(lv) > 100.0_dp .or. sigma_u <= 0.0_dp) then
                f = huge(1.0_dp)
                g = 0.0_dp
                return
            end if
            sig = exp(lv)
            z = regarch_obs(t) / max(sig, 1.0e-8_dp)
            z2 = z*z
            meas_hat = xi + phi*lv + tau1*z + tau2*(z2 - 1.0_dp)
            resid = regarch_log_range(t) - meas_hat
            u = resid / sigma_u
            f = f + log_sqrt_2pi + lv + 0.5_dp*z2 + log_sqrt_2pi + log(sigma_u) + 0.5_dp*u*u

            mcoeff = phi - tau1*z - 2.0_dp*tau2*z2
            factor_lv = 1.0_dp - z2 - resid*mcoeff/sigma_u**2
            grad = grad + factor_lv*d
            grad(5) = grad(5) - resid / sigma_u**2
            grad(6) = grad(6) - resid*lv / sigma_u**2
            grad(7) = grad(7) - resid*z / sigma_u**2
            grad(8) = grad(8) - resid*(z2 - 1.0_dp) / sigma_u**2
            grad(9) = grad(9) + 1.0_dp/sigma_u - resid**2/sigma_u**3

            dz = -z*d
            dm = mcoeff*d
            dm(5) = dm(5) + 1.0_dp
            dm(6) = dm(6) + lv
            dm(7) = dm(7) + z
            dm(8) = dm(8) + z2 - 1.0_dp
            du = -dm / sigma_u
            du(9) = du(9) - resid / sigma_u**2

            dnext = beta*d + alpha*du + gamma*dz
            dnext(1) = dnext(1) + 1.0_dp
            dnext(2) = dnext(2) + u
            dnext(3) = dnext(3) + z
            dnext(4) = dnext(4) + lv
            lv = omega + beta*lv + alpha*u + gamma*z
            d = dnext
        end do

        g = grad
        g(4) = g(4) * beta * (1.0_dp - beta)
        g(9) = g(9) * sigma_u
        f = f / real(regarch_nobs, dp)
        g = g / real(regarch_nobs, dp)
    end subroutine rgarch_meas_obj

    subroutine rgarch_meas_transform(p, omega, alpha, gamma, beta, xi, phi, tau1, tau2, sigma_u)
        real(dp), intent(in)  :: p(:)
        real(dp), intent(out) :: omega, alpha, gamma, beta, xi, phi, tau1, tau2, sigma_u

        omega = p(1)
        alpha = p(2)
        gamma = p(3)
        beta = 1.0_dp / (1.0_dp + exp(-p(4)))
        beta = min(max(beta, 1.0e-6_dp), 1.0_dp - 1.0e-6_dp)
        xi = p(5)
        phi = p(6)
        tau1 = p(7)
        tau2 = p(8)
        sigma_u = exp(p(9))
        sigma_u = max(sigma_u, 1.0e-6_dp)
    end subroutine rgarch_meas_transform

    subroutine rgarch_meas_inv_transform(omega, alpha, gamma, beta, xi, phi, tau1, tau2, sigma_u, p)
        real(dp), intent(in)  :: omega, alpha, gamma, beta, xi, phi, tau1, tau2, sigma_u
        real(dp), intent(out) :: p(rgarch_meas_np)
        real(dp) :: b

        b = min(max(beta, 1.0e-6_dp), 1.0_dp - 1.0e-6_dp)
        p(1) = omega
        p(2) = alpha
        p(3) = gamma
        p(4) = log(b / (1.0_dp - b))
        p(5) = xi
        p(6) = phi
        p(7) = tau1
        p(8) = tau2
        p(9) = log(max(sigma_u, 1.0e-6_dp))
    end subroutine rgarch_meas_inv_transform

    subroutine nagarch_range_set_data(y, range_var)
        real(dp), intent(in) :: y(:), range_var(:)

        nag_range_nobs = size(y)
        if (size(range_var) /= nag_range_nobs) then
            print '(A)', "garch_fit_mod: range_var length mismatch in fit_nagarch_range"
            error stop
        end if
        if (allocated(nag_range_obs)) deallocate(nag_range_obs)
        if (allocated(nag_range_x)) deallocate(nag_range_x)
        allocate(nag_range_obs(nag_range_nobs), nag_range_x(nag_range_nobs))
        nag_range_obs = y
        nag_range_x = range_var
    end subroutine nagarch_range_set_data

    subroutine nagarch_range_obj(p, np, f, g)
        integer,  intent(in)  :: np
        real(dp), intent(in)  :: p(np)
        real(dp), intent(out) :: f, g(np)

        real(dp) :: omega, alpha, gamma, beta, theta, h, sqrth, r, kappa
        real(dp) :: moment, dmoment, denom, numer, aa, factor
        real(dp) :: dh_dom, dh_dal, dh_dga, dh_dbe, dh_dth
        real(dp) :: grad_om, grad_al, grad_ga, grad_be, grad_th, xbar
        integer :: t

        call nagarch_range_transform(p, omega, alpha, gamma, beta, theta)
        call nagarch_shift_moment_local(theta, moment, dmoment)
        xbar = sum(nag_range_x) / real(nag_range_nobs, dp)
        denom = max(1.0_dp - alpha*moment - beta, 1.0e-8_dp)
        numer = omega + gamma*xbar
        h = max(numer / denom, sum(nag_range_obs**2)/real(nag_range_nobs, dp), 1.0e-12_dp)
        dh_dom = 1.0_dp / denom
        dh_dal = numer * moment / denom**2
        dh_dga = xbar / denom
        dh_dbe = numer / denom**2
        dh_dth = numer * alpha * dmoment / denom**2

        f = real(nag_range_nobs, dp) * log_sqrt_2pi
        grad_om = 0.0_dp
        grad_al = 0.0_dp
        grad_ga = 0.0_dp
        grad_be = 0.0_dp
        grad_th = 0.0_dp
        do t = 1, nag_range_nobs
            if (h <= 0.0_dp .or. h /= h .or. h > 1.0e100_dp) then
                f = huge(1.0_dp)
                g = 0.0_dp
                return
            end if
            h = max(h, 1.0e-12_dp)
            factor = 1.0_dp/h - nag_range_obs(t)**2/h**2
            f = f + 0.5_dp*(log(h) + nag_range_obs(t)**2/h)
            grad_om = grad_om + 0.5_dp*factor*dh_dom
            grad_al = grad_al + 0.5_dp*factor*dh_dal
            grad_ga = grad_ga + 0.5_dp*factor*dh_dga
            grad_be = grad_be + 0.5_dp*factor*dh_dbe
            grad_th = grad_th + 0.5_dp*factor*dh_dth
            sqrth = sqrt(h)
            r = nag_range_obs(t) - theta*sqrth
            kappa = beta - alpha*theta*r/sqrth
            dh_dom = 1.0_dp + kappa*dh_dom
            dh_dal = r**2 + kappa*dh_dal
            dh_dga = nag_range_x(t) + kappa*dh_dga
            dh_dbe = h + kappa*dh_dbe
            dh_dth = -2.0_dp*alpha*r*sqrth + kappa*dh_dth
            h = omega + alpha*r**2 + gamma*nag_range_x(t) + beta*h
        end do

        aa = alpha * moment
        g(1) = grad_om * omega
        g(2) = grad_al*alpha*(1.0_dp - aa) - grad_ga*gamma*aa - grad_be*beta*aa
        g(3) = -grad_al*alpha*gamma + grad_ga*gamma*(1.0_dp - gamma) - grad_be*beta*gamma
        g(4) = -grad_al*alpha*beta - grad_ga*gamma*beta + grad_be*beta*(1.0_dp - beta)
        g(5) = grad_al*(-alpha*dmoment/moment) + grad_th
        f = f / real(nag_range_nobs, dp)
        g = g / real(nag_range_nobs, dp)
    end subroutine nagarch_range_obj

    subroutine nagarch_range_transform(p, omega, alpha, gamma, beta, theta)
        real(dp), intent(in) :: p(:)
        real(dp), intent(out) :: omega, alpha, gamma, beta, theta
        real(dp) :: e2, e3, e4, s, aa, shift_moment, dmoment_unused

        omega = exp(p(1))
        theta = p(5)
        call nagarch_shift_moment_local(theta, shift_moment, dmoment_unused)
        e2 = exp(p(2))
        e3 = exp(p(3))
        e4 = exp(p(4))
        s = 1.0_dp + e2 + e3 + e4
        aa = e2 / s
        gamma = e3 / s
        beta = e4 / s
        alpha = aa / shift_moment
    end subroutine nagarch_range_transform

    subroutine nagarch_range_inv_transform(omega, alpha, gamma, beta, theta, p)
        real(dp), intent(in) :: omega, alpha, gamma, beta, theta
        real(dp), intent(out) :: p(nagarch_range_np)
        real(dp) :: aa, slack, shift_moment, dmoment_unused

        call nagarch_shift_moment_local(theta, shift_moment, dmoment_unused)
        aa = alpha * shift_moment
        slack = max(1.0_dp - aa - gamma - beta, 1.0e-8_dp)
        p(1) = log(max(omega, 1.0e-12_dp))
        p(2) = log(max(aa, 1.0e-12_dp) / slack)
        p(3) = log(max(gamma, 1.0e-12_dp) / slack)
        p(4) = log(max(beta, 1.0e-12_dp) / slack)
        p(5) = theta
    end subroutine nagarch_range_inv_transform

    real(dp) function nagarch_range_initial_variance(y, range_var, params)
        real(dp), intent(in) :: y(:), range_var(:)
        type(garch_params_t), intent(in) :: params
        nagarch_range_initial_variance = nagarch_range_initial_variance_raw(y, range_var, params%omega, &
            params%alpha, params%gamma, params%beta, params%theta)
    end function nagarch_range_initial_variance

    real(dp) function nagarch_range_initial_variance_raw(y, range_var, omega, alpha, gamma, beta, theta)
        real(dp), intent(in) :: y(:), range_var(:), omega, alpha, gamma, beta, theta
        real(dp) :: shift_moment, dmoment_unused, denom

        call nagarch_shift_moment_local(theta, shift_moment, dmoment_unused)
        denom = max(1.0_dp - alpha*shift_moment - beta, 1.0e-8_dp)
        nagarch_range_initial_variance_raw = max((omega + gamma*sum(range_var)/real(size(range_var), dp)) / denom, &
            sum(y**2)/real(size(y), dp), 1.0e-12_dp)
    end function nagarch_range_initial_variance_raw

    subroutine fgarch_twist_range_obj(p, np, f, g)
        integer,  intent(in)  :: np
        real(dp), intent(in)  :: p(np)
        real(dp), intent(out) :: f, g(np)

        real(dp) :: omega, alpha, gamma, beta, theta, twist, h, sqrth, z, u, sgn, q, q2
        real(dp) :: moment, dmth, dmc, denom, numer, aa, factor, xbar
        real(dp) :: dh_dom, dh_dal, dh_dga, dh_dbe, dh_dth, dh_dtw
        real(dp) :: dhn_dom, dhn_dal, dhn_dga, dhn_dbe, dhn_dth, dhn_dtw
        real(dp) :: dqn_dom, dqn_dal, dqn_dga, dqn_dbe, dqn_dth, dqn_dtw
        real(dp) :: grad_om, grad_al, grad_ga, grad_be, grad_th, grad_tw
        real(dp) :: h_old, z_over_h
        integer :: t

        call fgarch_twist_range_transform(p, omega, alpha, gamma, beta, theta, twist)
        call fgarch_twist_moment_local(theta, twist, moment, dmth, dmc)
        xbar = sum(nag_range_x) / real(nag_range_nobs, dp)
        denom = max(1.0_dp - alpha*moment - beta, 1.0e-8_dp)
        numer = omega + gamma*xbar
        h = max(numer / denom, sum(nag_range_obs**2)/real(nag_range_nobs, dp), 1.0e-12_dp)
        dh_dom = 1.0_dp / denom
        dh_dal = numer * moment / denom**2
        dh_dga = xbar / denom
        dh_dbe = numer / denom**2
        dh_dth = numer * alpha * dmth / denom**2
        dh_dtw = numer * alpha * dmc / denom**2

        f = real(nag_range_nobs, dp) * log_sqrt_2pi
        grad_om = 0.0_dp
        grad_al = 0.0_dp
        grad_ga = 0.0_dp
        grad_be = 0.0_dp
        grad_th = 0.0_dp
        grad_tw = 0.0_dp
        do t = 1, nag_range_nobs
            if (h <= 0.0_dp .or. h /= h .or. h > 1.0e100_dp) then
                f = huge(1.0_dp)
                g = 0.0_dp
                return
            end if
            h = max(h, 1.0e-12_dp)
            factor = 1.0_dp/h - nag_range_obs(t)**2/h**2
            f = f + 0.5_dp*(log(h) + nag_range_obs(t)**2/h)
            grad_om = grad_om + 0.5_dp*factor*dh_dom
            grad_al = grad_al + 0.5_dp*factor*dh_dal
            grad_ga = grad_ga + 0.5_dp*factor*dh_dga
            grad_be = grad_be + 0.5_dp*factor*dh_dbe
            grad_th = grad_th + 0.5_dp*factor*dh_dth
            grad_tw = grad_tw + 0.5_dp*factor*dh_dtw

            h_old = h
            sqrth = sqrt(h_old)
            z = nag_range_obs(t) / sqrth
            u = z - theta
            sgn = merge(1.0_dp, -1.0_dp, u >= 0.0_dp)
            q = abs(u) - twist*u
            q2 = q*q
            z_over_h = z / h_old

            dqn_dom = 2.0_dp*q*(sgn - twist)*(-0.5_dp*z_over_h*dh_dom)
            dqn_dal = 2.0_dp*q*(sgn - twist)*(-0.5_dp*z_over_h*dh_dal)
            dqn_dga = 2.0_dp*q*(sgn - twist)*(-0.5_dp*z_over_h*dh_dga)
            dqn_dbe = 2.0_dp*q*(sgn - twist)*(-0.5_dp*z_over_h*dh_dbe)
            dqn_dth = 2.0_dp*q*((sgn - twist)*(-0.5_dp*z_over_h*dh_dth - 1.0_dp))
            dqn_dtw = 2.0_dp*q*((sgn - twist)*(-0.5_dp*z_over_h*dh_dtw) - u)

            dhn_dom = 1.0_dp + alpha*(q2*dh_dom + h_old*dqn_dom) + beta*dh_dom
            dhn_dal = h_old*q2 + alpha*(q2*dh_dal + h_old*dqn_dal) + beta*dh_dal
            dhn_dga = nag_range_x(t) + alpha*(q2*dh_dga + h_old*dqn_dga) + beta*dh_dga
            dhn_dbe = h_old + alpha*(q2*dh_dbe + h_old*dqn_dbe) + beta*dh_dbe
            dhn_dth = alpha*(q2*dh_dth + h_old*dqn_dth) + beta*dh_dth
            dhn_dtw = alpha*(q2*dh_dtw + h_old*dqn_dtw) + beta*dh_dtw
            h = omega + alpha*h_old*q2 + gamma*nag_range_x(t) + beta*h_old
            dh_dom = dhn_dom
            dh_dal = dhn_dal
            dh_dga = dhn_dga
            dh_dbe = dhn_dbe
            dh_dth = dhn_dth
            dh_dtw = dhn_dtw
        end do

        aa = alpha * moment
        g(1) = grad_om * omega
        g(2) = grad_al*alpha*(1.0_dp - aa) - grad_ga*gamma*aa - grad_be*beta*aa
        g(3) = -grad_al*alpha*gamma + grad_ga*gamma*(1.0_dp - gamma) - grad_be*beta*gamma
        g(4) = -grad_al*alpha*beta - grad_ga*gamma*beta + grad_be*beta*(1.0_dp - beta)
        g(5) = grad_al*(-alpha*dmth/moment) + grad_th
        g(6) = grad_al*(-alpha*dmc/moment) + grad_tw
        f = f / real(nag_range_nobs, dp)
        g = g / real(nag_range_nobs, dp)
    end subroutine fgarch_twist_range_obj

    subroutine fgarch_twist_range_transform(p, omega, alpha, gamma, beta, theta, twist)
        real(dp), intent(in) :: p(:)
        real(dp), intent(out) :: omega, alpha, gamma, beta, theta, twist
        real(dp) :: e2, e3, e4, s, aa, moment, dmth, dmc

        omega = exp(p(1))
        theta = p(5)
        twist = p(6)
        call fgarch_twist_moment_local(theta, twist, moment, dmth, dmc)
        e2 = exp(p(2))
        e3 = exp(p(3))
        e4 = exp(p(4))
        s = 1.0_dp + e2 + e3 + e4
        aa = e2 / s
        gamma = e3 / s
        beta = e4 / s
        alpha = aa / moment
    end subroutine fgarch_twist_range_transform

    subroutine fgarch_twist_range_inv_transform(omega, alpha, gamma, beta, theta, twist, p)
        real(dp), intent(in) :: omega, alpha, gamma, beta, theta, twist
        real(dp), intent(out) :: p(fgarch_twist_range_np)
        real(dp) :: aa, slack, moment, dmth, dmc

        call fgarch_twist_moment_local(theta, twist, moment, dmth, dmc)
        aa = alpha * moment
        slack = max(1.0_dp - aa - gamma - beta, 1.0e-8_dp)
        p(1) = log(max(omega, 1.0e-12_dp))
        p(2) = log(max(aa, 1.0e-12_dp) / slack)
        p(3) = log(max(gamma, 1.0e-12_dp) / slack)
        p(4) = log(max(beta, 1.0e-12_dp) / slack)
        p(5) = theta
        p(6) = twist
    end subroutine fgarch_twist_range_inv_transform

    real(dp) function fgarch_twist_range_initial_variance(y, range_var, params)
        real(dp), intent(in) :: y(:), range_var(:)
        type(garch_params_t), intent(in) :: params
        real(dp) :: moment, dmth, dmc, denom

        call fgarch_twist_moment_local(params%theta, params%twist, moment, dmth, dmc)
        denom = max(1.0_dp - params%alpha*moment - params%beta, 1.0e-8_dp)
        fgarch_twist_range_initial_variance = max((params%omega + params%gamma*sum(range_var)/real(size(range_var), dp)) / &
            denom, sum(y**2)/real(size(y), dp), 1.0e-12_dp)
    end function fgarch_twist_range_initial_variance

    subroutine fgarch_twist_moment_local(theta, twist, moment, dmth, dmc)
        real(dp), intent(in) :: theta, twist
        real(dp), intent(out) :: moment, dmth, dmc
        real(dp) :: Phi, ph, A, B, dA, dB, one_minus, one_plus

        Phi = 0.5_dp * (1.0_dp + erf(theta / sqrt2))
        ph  = exp(-0.5_dp * theta**2 - log_sqrt_2pi)
        A = (1.0_dp + theta**2) * (1.0_dp - Phi) - theta * ph
        B = (1.0_dp + theta**2) * Phi + theta * ph
        dA = 2.0_dp * (theta*(1.0_dp - Phi) - ph)
        dB = 2.0_dp * (theta*Phi + ph)
        one_minus = 1.0_dp - twist
        one_plus  = 1.0_dp + twist
        moment = one_minus**2 * A + one_plus**2 * B
        dmth = one_minus**2 * dA + one_plus**2 * dB
        dmc = -2.0_dp*one_minus*A + 2.0_dp*one_plus*B
    end subroutine fgarch_twist_moment_local

    subroutine nagarch_shift_moment_local(theta, moment, dmoment)
        real(dp), intent(in) :: theta
        real(dp), intent(out) :: moment, dmoment

        moment = 1.0_dp + theta**2
        dmoment = 2.0_dp*theta
    end subroutine nagarch_shift_moment_local

    subroutine aewma_nag_obj(p, np, f, g)
        integer,  intent(in)  :: np
        real(dp), intent(in)  :: p(np)
        real(dp), intent(out) :: f, g(np)

        call aewma_nag_value_grad(p, f, g)
    end subroutine aewma_nag_obj

    subroutine aewma_twist_obj(p, np, f, g)
        integer,  intent(in)  :: np
        real(dp), intent(in)  :: p(np)
        real(dp), intent(out) :: f, g(np)

        call aewma_twist_value_grad(p, f, g)
    end subroutine aewma_twist_obj

    subroutine aewma_nag_value_grad(p, f, g)
        real(dp), intent(in)  :: p(:)
        real(dp), intent(out) :: f, g(:)
        real(dp) :: lambda, theta, scale, h, z, u, nic, mult, fsum, yt2, dfdh
        real(dp) :: dh(aewma_nag_np), dnic(aewma_nag_np), dmult(aewma_nag_np)
        real(dp) :: dlambda(aewma_nag_np), dtheta(aewma_nag_np), dscale(aewma_nag_np)
        integer :: t

        call aewma_nag_transform(p, lambda, theta, scale)
        dlambda = 0.0_dp
        dtheta = 0.0_dp
        dscale = 0.0_dp
        dlambda(1) = lambda * (1.0_dp - lambda)
        dtheta(2) = 1.0_dp
        dscale(3) = scale
        h = sample_variance(ewma_obs)
        dh = 0.0_dp
        fsum = 0.0_dp
        g = 0.0_dp

        do t = 1, ewma_nobs
            if (h <= 0.0_dp .or. h /= h .or. h > 1.0e100_dp) then
                f = huge(1.0_dp)
                g = 0.0_dp
                return
            end if
            h = max(h, 1.0e-12_dp)
            yt2 = ewma_obs(t)**2
            z = ewma_obs(t) / sqrt(h)
            fsum = fsum + log_sqrt_2pi + 0.5_dp*(log(h) + yt2/h)
            dfdh = 0.5_dp*(1.0_dp - yt2/h) / h
            g = g + dfdh*dh

            u = z - theta
            nic = u**2
            mult = lambda + (1.0_dp - lambda)*scale*nic
            dnic = 2.0_dp*u*(-0.5_dp*z*dh/h - dtheta)
            dmult = dlambda*(1.0_dp - scale*nic) + &
                (1.0_dp - lambda)*(dscale*nic + scale*dnic)
            dh = dh*mult + h*dmult
            h = h*mult
        end do

        f = fsum / real(ewma_nobs, dp)
        g = g / real(ewma_nobs, dp)
    end subroutine aewma_nag_value_grad

    subroutine aewma_twist_value_grad(p, f, g)
        real(dp), intent(in)  :: p(:)
        real(dp), intent(out) :: f, g(:)
        real(dp) :: lambda, theta, twist, scale, h, z, u, sgn, q, nic, mult, fsum, yt2, dfdh
        real(dp) :: dh(aewma_twist_np), du(aewma_twist_np), dq(aewma_twist_np), dnic(aewma_twist_np)
        real(dp) :: dmult(aewma_twist_np)
        real(dp) :: dlambda(aewma_twist_np), dtheta(aewma_twist_np), dtwist(aewma_twist_np), dscale(aewma_twist_np)
        integer :: t

        call aewma_twist_transform(p, lambda, theta, twist, scale)
        dlambda = 0.0_dp
        dtheta = 0.0_dp
        dtwist = 0.0_dp
        dscale = 0.0_dp
        dlambda(1) = lambda * (1.0_dp - lambda)
        dtheta(2) = 1.0_dp
        dtwist(3) = 1.0_dp
        dscale(4) = scale
        h = sample_variance(ewma_obs)
        dh = 0.0_dp
        fsum = 0.0_dp
        g = 0.0_dp

        do t = 1, ewma_nobs
            if (h <= 0.0_dp .or. h /= h .or. h > 1.0e100_dp) then
                f = huge(1.0_dp)
                g = 0.0_dp
                return
            end if
            h = max(h, 1.0e-12_dp)
            yt2 = ewma_obs(t)**2
            z = ewma_obs(t) / sqrt(h)
            fsum = fsum + log_sqrt_2pi + 0.5_dp*(log(h) + yt2/h)
            dfdh = 0.5_dp*(1.0_dp - yt2/h) / h
            g = g + dfdh*dh

            u = z - theta
            sgn = merge(1.0_dp, -1.0_dp, u >= 0.0_dp)
            q = abs(u) - twist*u
            nic = q**2
            mult = lambda + (1.0_dp - lambda)*scale*nic
            du = -0.5_dp*z*dh/h - dtheta
            dq = (sgn - twist)*du - dtwist*u
            dnic = 2.0_dp*q*dq
            dmult = dlambda*(1.0_dp - scale*nic) + &
                (1.0_dp - lambda)*(dscale*nic + scale*dnic)
            dh = dh*mult + h*dmult
            h = h*mult
        end do

        f = fsum / real(ewma_nobs, dp)
        g = g / real(ewma_nobs, dp)
    end subroutine aewma_twist_value_grad

    subroutine aewma_nag_transform(p, lambda, theta, scale)
        real(dp), intent(in)  :: p(:)
        real(dp), intent(out) :: lambda, theta, scale

        lambda = 1.0_dp / (1.0_dp + exp(-p(1)))
        lambda = min(max(lambda, 1.0e-6_dp), 1.0_dp - 1.0e-6_dp)
        theta = p(2)
        scale = exp(p(3))
    end subroutine aewma_nag_transform

    subroutine aewma_twist_transform(p, lambda, theta, twist, scale)
        real(dp), intent(in)  :: p(:)
        real(dp), intent(out) :: lambda, theta, twist, scale

        lambda = 1.0_dp / (1.0_dp + exp(-p(1)))
        lambda = min(max(lambda, 1.0e-6_dp), 1.0_dp - 1.0e-6_dp)
        theta = p(2)
        twist = p(3)
        scale = exp(p(4))
    end subroutine aewma_twist_transform

    subroutine aewma_nag_inv_transform(lambda, theta, scale, p)
        real(dp), intent(in)  :: lambda, theta, scale
        real(dp), intent(out) :: p(aewma_nag_np)
        real(dp) :: lam

        lam = min(max(lambda, 1.0e-6_dp), 1.0_dp - 1.0e-6_dp)
        p(1) = log(lam / (1.0_dp - lam))
        p(2) = theta
        p(3) = log(max(scale, 1.0e-12_dp))
    end subroutine aewma_nag_inv_transform

    subroutine aewma_twist_inv_transform(lambda, theta, twist, scale, p)
        real(dp), intent(in)  :: lambda, theta, twist, scale
        real(dp), intent(out) :: p(aewma_twist_np)
        real(dp) :: lam

        lam = min(max(lambda, 1.0e-6_dp), 1.0_dp - 1.0e-6_dp)
        p(1) = log(lam / (1.0_dp - lam))
        p(2) = theta
        p(3) = twist
        p(4) = log(max(scale, 1.0e-12_dp))
    end subroutine aewma_twist_inv_transform

    real(dp) function sample_variance(y)
        real(dp), intent(in) :: y(:)
        sample_variance = max(sum(y**2) / real(size(y), dp), 1.0e-12_dp)
    end function sample_variance

    subroutine moments(zz, n, skew, ekurt)
        integer,  intent(in)  :: n
        real(dp), intent(in)  :: zz(n)
        real(dp), intent(out) :: skew, ekurt
        real(dp) :: rn, zm, zv, zs, zk, dz
        integer :: t

        rn = real(n, dp)
        zm = sum(zz) / rn
        zv = 0.0_dp
        zs = 0.0_dp
        zk = 0.0_dp
        do t = 1, n
            dz = zz(t) - zm
            zv = zv + dz**2
            zs = zs + dz**3
            zk = zk + dz**4
        end do
        zv = zv / rn
        skew = (zs / rn) / zv**1.5_dp
        ekurt = (zk / rn) / zv**2 - 3.0_dp
    end subroutine moments

end module garch_fit_mod

! Daily volatility forecasts built from realized volatility measures.
!
! The first model implemented here is an EWMA-affine variance forecast:
!     m_t = lambda*m_{t-1} + (1-lambda)*x_{t-1}
!     h_t = a + b*m_t
! where x_t is a daily realized variance measure computed from intraday data.

module realized_vol_forecast_mod
    use kind_mod, only: dp
    use math_const_mod, only: log_sqrt_2pi
    use bfgs_mod, only: bfgs_minimize
    implicit none
    private

    real(dp), parameter :: min_var = 1.0e-12_dp

    type, public :: ewma_affine_result_t
        real(dp) :: lambda = 0.94_dp
        real(dp) :: a = 0.0_dp
        real(dp) :: b = 1.0_dp
        real(dp) :: loglik = -huge(1.0_dp)
        integer :: niter = 0
        logical :: converged = .false.
    end type ewma_affine_result_t

    type, public :: affine_variance_result_t
        real(dp) :: a = 0.0_dp
        real(dp) :: b = 1.0_dp
        real(dp) :: loglik = -huge(1.0_dp)
        integer :: niter = 0
        logical :: converged = .false.
    end type affine_variance_result_t

    type, public :: affine2_variance_result_t
        real(dp) :: a = 0.0_dp
        real(dp) :: b1 = 1.0_dp
        real(dp) :: b2 = 1.0_dp
        real(dp) :: loglik = -huge(1.0_dp)
        integer :: niter = 0
        logical :: converged = .false.
    end type affine2_variance_result_t

    type, public :: har_variance_result_t
        real(dp) :: coef(4) = 0.0_dp
        real(dp) :: loglik = -huge(1.0_dp)
        integer :: niter = 0
        logical :: converged = .false.
    end type har_variance_result_t

    type, public :: harx_variance_result_t
        real(dp) :: coef(5) = 0.0_dp
        real(dp) :: loglik = -huge(1.0_dp)
        integer :: niter = 0
        logical :: converged = .false.
    end type harx_variance_result_t

    type, public :: harx_lev_variance_result_t
        real(dp) :: coef(8) = 0.0_dp
        real(dp) :: loglik = -huge(1.0_dp)
        integer :: niter = 0
        logical :: converged = .false.
    end type harx_lev_variance_result_t

    type, public :: log_har_variance_result_t
        real(dp) :: coef(4) = 0.0_dp
        real(dp) :: loglik = -huge(1.0_dp)
        integer :: niter = 0
        logical :: converged = .false.
    end type log_har_variance_result_t

    type, public :: sqrt_har_variance_result_t
        real(dp) :: coef(4) = 0.0_dp
        real(dp) :: loglik = -huge(1.0_dp)
        integer :: niter = 0
        logical :: converged = .false.
    end type sqrt_har_variance_result_t

    type, public :: harq_variance_result_t
        real(dp) :: coef(5) = 0.0_dp
        real(dp) :: loglik = -huge(1.0_dp)
        integer :: niter = 0
        logical :: converged = .false.
    end type harq_variance_result_t

    type, public :: harj_variance_result_t
        real(dp) :: coef(7) = 0.0_dp
        real(dp) :: loglik = -huge(1.0_dp)
        integer :: niter = 0
        logical :: converged = .false.
    end type harj_variance_result_t

    type, public :: harqj_variance_result_t
        real(dp) :: coef(8) = 0.0_dp
        real(dp) :: loglik = -huge(1.0_dp)
        integer :: niter = 0
        logical :: converged = .false.
    end type harqj_variance_result_t

    type, public :: semivar_har_variance_result_t
        real(dp) :: coef(7) = 0.0_dp
        real(dp) :: loglik = -huge(1.0_dp)
        integer :: niter = 0
        logical :: converged = .false.
    end type semivar_har_variance_result_t

    type, public :: midas_variance_result_t
        real(dp) :: a = 0.0_dp
        real(dp) :: b = 1.0_dp
        real(dp) :: theta1 = 1.0_dp
        real(dp) :: theta2 = 2.0_dp
        real(dp) :: loglik = -huge(1.0_dp)
        integer :: niter = 0
        logical :: converged = .false.
    end type midas_variance_result_t

    type, public :: har_negret_variance_result_t
        real(dp) :: coef(5) = 0.0_dp
        real(dp) :: loglik = -huge(1.0_dp)
        integer :: niter = 0
        logical :: converged = .false.
    end type har_negret_variance_result_t

    type, public :: har_lev_variance_result_t
        real(dp) :: coef(7) = 0.0_dp
        real(dp) :: loglik = -huge(1.0_dp)
        integer :: niter = 0
        logical :: converged = .false.
    end type har_lev_variance_result_t

    type, public :: heavy_variance_result_t
        real(dp) :: omega = 0.0_dp
        real(dp) :: alpha = 0.0_dp
        real(dp) :: beta = 0.0_dp
        real(dp) :: omega_rm = 0.0_dp
        real(dp) :: alpha_rm = 0.0_dp
        real(dp) :: beta_rm = 0.0_dp
        real(dp) :: loglik = -huge(1.0_dp)
        real(dp) :: loglik_rm = -huge(1.0_dp)
        integer :: niter = 0
        logical :: converged = .false.
    end type heavy_variance_result_t

    real(dp), allocatable, save :: obj_y(:), obj_x(:)
    real(dp), allocatable, save :: obj_x2(:), obj_x3(:)
    integer, save :: obj_ntrain = 0
    integer, save :: obj_k_lag = 22
    logical, save :: obj_fit_lambda = .false.
    real(dp), save :: obj_fixed_lambda = 0.94_dp

    public :: fit_ewma_affine_variance
    public :: fit_affine_variance
    public :: fit_affine2_variance
    public :: fit_har_variance
    public :: fit_harx_variance
    public :: fit_harx_lev_variance
    public :: fit_log_har_variance
    public :: fit_sqrt_har_variance
    public :: fit_harq_variance
    public :: fit_harj_variance
    public :: fit_harqj_variance
    public :: fit_har_negret_variance
    public :: fit_har_lev_variance
    public :: fit_semivar_har_variance
    public :: fit_midas_variance
    public :: fit_heavy_variance
    public :: har_variance_path
    public :: harx_variance_path
    public :: harx_lev_variance_path
    public :: log_har_variance_path
    public :: sqrt_har_variance_path
    public :: harq_variance_path
    public :: harj_variance_path
    public :: harqj_variance_path
    public :: har_negret_variance_path
    public :: har_lev_variance_path
    public :: semivar_har_variance_path
    public :: midas_variance_path
    public :: heavy_variance_path
    public :: affine_variance_path
    public :: affine2_variance_path
    public :: ewma_affine_variance_path
    public :: gaussian_variance_loglik
    public :: qlike_loss

contains

    subroutine fit_ewma_affine_variance(y, x, ntrain, fit_lambda, lambda_fixed, max_iter, gtol, result, h)
        ! Fit h_t = a + b*EWMA(x) by Gaussian likelihood on y(1:ntrain).
        real(dp), intent(in) :: y(:), x(:), lambda_fixed, gtol
        integer, intent(in) :: ntrain, max_iter
        logical, intent(in) :: fit_lambda
        type(ewma_affine_result_t), intent(out) :: result
        real(dp), intent(out) :: h(:)
        real(dp), allocatable :: p(:), grad(:)
        real(dp) :: fopt, mean_y2, mean_x
        integer :: np, niter
        logical :: converged

        if (size(y) /= size(x) .or. size(h) /= size(y)) error stop "fit_ewma_affine_variance: array sizes differ"
        if (ntrain < 5 .or. ntrain > size(y)) error stop "fit_ewma_affine_variance: invalid ntrain"
        if (max_iter < 1) error stop "fit_ewma_affine_variance: max_iter must be positive"

        if (allocated(obj_y)) deallocate(obj_y)
        if (allocated(obj_x)) deallocate(obj_x)
        allocate(obj_y(size(y)), obj_x(size(x)))
        obj_y = y
        obj_x = max(x, min_var)
        obj_ntrain = ntrain
        obj_fit_lambda = fit_lambda
        obj_fixed_lambda = min(max(lambda_fixed, 0.0001_dp), 0.9999_dp)

        np = merge(3, 2, fit_lambda)
        allocate(p(np), grad(np))
        mean_y2 = max(sum(y(1:ntrain)**2) / real(ntrain, dp), min_var)
        mean_x = max(sum(obj_x(1:ntrain)) / real(ntrain, dp), min_var)
        p(1) = log(0.1_dp*mean_y2)
        p(2) = log(max(0.9_dp*mean_y2 / mean_x, min_var))
        if (fit_lambda) p(3) = lambda_inverse(0.94_dp)

        call bfgs_minimize(ewma_affine_obj, p, np, max_iter, gtol, fopt, niter, converged)
        result%a = exp(p(1))
        result%b = exp(p(2))
        if (fit_lambda) then
            result%lambda = lambda_transform(p(3))
        else
            result%lambda = obj_fixed_lambda
        end if
        call ewma_affine_variance_path(x, result%lambda, result%a, result%b, ntrain, h)
        result%loglik = gaussian_variance_loglik(y(1:ntrain), h(1:ntrain))
        result%niter = niter
        result%converged = converged
        deallocate(p, grad)
    end subroutine fit_ewma_affine_variance

    subroutine fit_affine_variance(y, x, ntrain, max_iter, gtol, result, h)
        ! Fit h_t = a + b*x_t by Gaussian likelihood with positive a and b.
        real(dp), intent(in) :: y(:), x(:), gtol
        integer, intent(in) :: ntrain, max_iter
        type(affine_variance_result_t), intent(out) :: result
        real(dp), intent(out) :: h(:)
        real(dp) :: p(2), fopt, mean_y2, mean_x
        integer :: niter
        logical :: converged

        if (size(y) /= size(x) .or. size(h) /= size(y)) error stop "fit_affine_variance: array sizes differ"
        if (ntrain < 5 .or. ntrain > size(y)) error stop "fit_affine_variance: invalid ntrain"
        if (max_iter < 1) error stop "fit_affine_variance: max_iter must be positive"
        if (allocated(obj_y)) deallocate(obj_y)
        if (allocated(obj_x)) deallocate(obj_x)
        allocate(obj_y(size(y)), obj_x(size(x)))
        obj_y = y
        obj_x = max(x, min_var)
        obj_ntrain = ntrain

        mean_y2 = max(sum(y(1:ntrain)**2) / real(ntrain, dp), min_var)
        mean_x = max(sum(obj_x(1:ntrain)) / real(ntrain, dp), min_var)
        p(1) = log(0.1_dp*mean_y2)
        p(2) = log(max(0.9_dp*mean_y2 / mean_x, min_var))
        call bfgs_minimize(affine_obj, p, 2, max_iter, gtol, fopt, niter, converged)
        result%a = exp(p(1))
        result%b = exp(p(2))
        call affine_variance_path(x, result%a, result%b, h)
        result%loglik = gaussian_variance_loglik(y(1:ntrain), h(1:ntrain))
        result%niter = niter
        result%converged = converged
    end subroutine fit_affine_variance

    subroutine fit_affine2_variance(y, x1, x2, ntrain, max_iter, gtol, result, h)
        ! Fit h_t = a + b1*x1_t + b2*x2_t by Gaussian likelihood with positive coefficients.
        real(dp), intent(in) :: y(:), x1(:), x2(:), gtol
        integer, intent(in) :: ntrain, max_iter
        type(affine2_variance_result_t), intent(out) :: result
        real(dp), intent(out) :: h(:)
        real(dp) :: p(3), fopt, mean_y2, mean_x1, mean_x2
        integer :: niter
        logical :: converged

        if (size(y) /= size(x1) .or. size(y) /= size(x2) .or. size(h) /= size(y)) then
            error stop "fit_affine2_variance: array sizes differ"
        end if
        if (ntrain < 5 .or. ntrain > size(y)) error stop "fit_affine2_variance: invalid ntrain"
        if (max_iter < 1) error stop "fit_affine2_variance: max_iter must be positive"
        if (allocated(obj_y)) deallocate(obj_y)
        if (allocated(obj_x)) deallocate(obj_x)
        if (allocated(obj_x2)) deallocate(obj_x2)
        allocate(obj_y(size(y)), obj_x(size(x1)), obj_x2(size(x2)))
        obj_y = y
        obj_x = max(x1, min_var)
        obj_x2 = max(x2, min_var)
        obj_ntrain = ntrain

        mean_y2 = max(sum(y(1:ntrain)**2) / real(ntrain, dp), min_var)
        mean_x1 = max(sum(obj_x(1:ntrain)) / real(ntrain, dp), min_var)
        mean_x2 = max(sum(obj_x2(1:ntrain)) / real(ntrain, dp), min_var)
        p(1) = log(0.1_dp*mean_y2)
        p(2) = log(max(0.45_dp*mean_y2 / mean_x1, min_var))
        p(3) = log(max(0.45_dp*mean_y2 / mean_x2, min_var))
        call bfgs_minimize(affine2_obj, p, 3, max_iter, gtol, fopt, niter, converged)
        result%a = exp(p(1))
        result%b1 = exp(p(2))
        result%b2 = exp(p(3))
        call affine2_variance_path(x1, x2, result%a, result%b1, result%b2, h)
        result%loglik = gaussian_variance_loglik(y(1:ntrain), h(1:ntrain))
        result%niter = niter
        result%converged = converged
    end subroutine fit_affine2_variance

    subroutine fit_har_variance(y, x, ntrain, result, h)
        ! Fit h_t = c + b1*x_{t-1} + b5*avg_5(x) + b22*avg_22(x) by OLS on y_t^2.
        real(dp), intent(in) :: y(:), x(:)
        integer, intent(in) :: ntrain
        type(har_variance_result_t), intent(out) :: result
        real(dp), intent(out) :: h(:)
        real(dp) :: p(4), fopt, mean_y2, mean_x
        integer :: niter
        logical :: converged

        if (size(y) /= size(x) .or. size(h) /= size(y)) error stop "fit_har_variance: array sizes differ"
        if (ntrain < 30 .or. ntrain > size(y)) error stop "fit_har_variance: invalid ntrain"
        if (allocated(obj_y)) deallocate(obj_y)
        if (allocated(obj_x)) deallocate(obj_x)
        allocate(obj_y(size(y)), obj_x(size(x)))
        obj_y = y
        obj_x = max(x, min_var)
        obj_ntrain = ntrain

        mean_y2 = max(sum(y(1:ntrain)**2) / real(ntrain, dp), min_var)
        mean_x = max(sum(obj_x(1:ntrain)) / real(ntrain, dp), min_var)
        p(1) = log(0.25_dp*mean_y2)
        p(2:4) = log(max(0.25_dp*mean_y2 / mean_x, min_var))
        call bfgs_minimize(har_obj, p, 4, 300, 1.0e-5_dp, fopt, niter, converged)
        result%coef = exp(p)
        result%niter = niter
        result%converged = converged
        call har_variance_path(x, result%coef, h)
        result%loglik = gaussian_variance_loglik(y(23:ntrain), h(23:ntrain))
    end subroutine fit_har_variance

    subroutine fit_harx_variance(y, x, x_exog, ntrain, result, h)
        ! Fit positive-coefficient HAR plus a nonnegative exogenous variance predictor.
        real(dp), intent(in) :: y(:), x(:), x_exog(:)
        integer, intent(in) :: ntrain
        type(harx_variance_result_t), intent(out) :: result
        real(dp), intent(out) :: h(:)
        real(dp) :: p(5), fopt, mean_y2, mean_x, mean_exog
        integer :: niter
        logical :: converged

        if (size(y) /= size(x) .or. size(y) /= size(x_exog) .or. size(h) /= size(y)) then
            error stop "fit_harx_variance: array sizes differ"
        end if
        if (ntrain < 30 .or. ntrain > size(y)) error stop "fit_harx_variance: invalid ntrain"
        if (allocated(obj_y)) deallocate(obj_y)
        if (allocated(obj_x)) deallocate(obj_x)
        if (allocated(obj_x2)) deallocate(obj_x2)
        allocate(obj_y(size(y)), obj_x(size(x)), obj_x2(size(x_exog)))
        obj_y = y
        obj_x = max(x, min_var)
        obj_x2 = max(x_exog, min_var)
        obj_ntrain = ntrain

        mean_y2 = max(sum(y(1:ntrain)**2) / real(ntrain, dp), min_var)
        mean_x = max(sum(obj_x(1:ntrain)) / real(ntrain, dp), min_var)
        mean_exog = max(sum(obj_x2(1:ntrain)) / real(ntrain, dp), min_var)
        p(1) = log(0.20_dp*mean_y2)
        p(2:4) = log(max(0.20_dp*mean_y2 / mean_x, min_var))
        p(5) = log(max(0.20_dp*mean_y2 / mean_exog, min_var))
        call bfgs_minimize(harx_obj, p, 5, 300, 1.0e-5_dp, fopt, niter, converged)
        result%coef = exp(p)
        result%niter = niter
        result%converged = converged
        call harx_variance_path(x, x_exog, result%coef, h)
        result%loglik = gaussian_variance_loglik(y(23:ntrain), h(23:ntrain))
    end subroutine fit_harx_variance

    subroutine fit_harx_lev_variance(y, x, x_exog, ntrain, result, h)
        ! Fit HAR plus exogenous variance and daily/weekly/monthly negative-return-squared terms.
        real(dp), intent(in) :: y(:), x(:), x_exog(:)
        integer, intent(in) :: ntrain
        type(harx_lev_variance_result_t), intent(out) :: result
        real(dp), intent(out) :: h(:)
        real(dp) :: p(8), fopt, mean_y2, mean_x, mean_exog, mean_neg
        integer :: niter
        logical :: converged

        if (size(y) /= size(x) .or. size(y) /= size(x_exog) .or. size(h) /= size(y)) then
            error stop "fit_harx_lev_variance: array sizes differ"
        end if
        if (ntrain < 30 .or. ntrain > size(y)) error stop "fit_harx_lev_variance: invalid ntrain"
        if (allocated(obj_y)) deallocate(obj_y)
        if (allocated(obj_x)) deallocate(obj_x)
        if (allocated(obj_x2)) deallocate(obj_x2)
        allocate(obj_y(size(y)), obj_x(size(x)), obj_x2(size(x_exog)))
        obj_y = y
        obj_x = max(x, min_var)
        obj_x2 = max(x_exog, min_var)
        obj_ntrain = ntrain

        mean_y2 = max(sum(y(1:ntrain)**2) / real(ntrain, dp), min_var)
        mean_x = max(sum(obj_x(1:ntrain)) / real(ntrain, dp), min_var)
        mean_exog = max(sum(obj_x2(1:ntrain)) / real(ntrain, dp), min_var)
        mean_neg = max(sum(min(y(1:ntrain), 0.0_dp)**2) / real(ntrain, dp), min_var)
        p(1) = log(0.125_dp*mean_y2)
        p(2:4) = log(max(0.125_dp*mean_y2 / mean_x, min_var))
        p(5) = log(max(0.125_dp*mean_y2 / mean_exog, min_var))
        p(6:8) = log(max(0.125_dp*mean_y2 / mean_neg, min_var))
        call bfgs_minimize(harx_lev_obj, p, 8, 300, 1.0e-5_dp, fopt, niter, converged)
        result%coef = exp(p)
        result%niter = niter
        result%converged = converged
        call harx_lev_variance_path(y, x, x_exog, result%coef, h)
        result%loglik = gaussian_variance_loglik(y(23:ntrain), h(23:ntrain))
    end subroutine fit_harx_lev_variance

    subroutine fit_log_har_variance(y, x, ntrain, result, h)
        ! Fit log h_t = c + b1*log(x_{t-1}) + b5*log(avg_5(x)) + b22*log(avg_22(x)).
        real(dp), intent(in) :: y(:), x(:)
        integer, intent(in) :: ntrain
        type(log_har_variance_result_t), intent(out) :: result
        real(dp), intent(out) :: h(:)
        real(dp) :: p(4), fopt, mean_y2
        integer :: niter
        logical :: converged

        if (size(y) /= size(x) .or. size(h) /= size(y)) error stop "fit_log_har_variance: array sizes differ"
        if (ntrain < 30 .or. ntrain > size(y)) error stop "fit_log_har_variance: invalid ntrain"
        if (allocated(obj_y)) deallocate(obj_y)
        if (allocated(obj_x)) deallocate(obj_x)
        allocate(obj_y(size(y)), obj_x(size(x)))
        obj_y = y
        obj_x = max(x, min_var)
        obj_ntrain = ntrain

        mean_y2 = max(sum(y(1:ntrain)**2) / real(ntrain, dp), min_var)
        p = 0.0_dp
        p(1) = log(mean_y2)
        call bfgs_minimize(log_har_obj, p, 4, 300, 1.0e-5_dp, fopt, niter, converged)
        result%coef = p
        result%niter = niter
        result%converged = converged
        call log_har_variance_path(x, result%coef, h)
        result%loglik = gaussian_variance_loglik(y(23:ntrain), h(23:ntrain))
    end subroutine fit_log_har_variance

    subroutine fit_sqrt_har_variance(y, x, ntrain, result, h)
        ! Fit sqrt(h_t) = c + b1*sqrt(x_{t-1}) + b5*sqrt(avg_5(x)) + b22*sqrt(avg_22(x)).
        real(dp), intent(in) :: y(:), x(:)
        integer, intent(in) :: ntrain
        type(sqrt_har_variance_result_t), intent(out) :: result
        real(dp), intent(out) :: h(:)
        real(dp) :: p(4), fopt, mean_y2
        integer :: niter
        logical :: converged

        if (size(y) /= size(x) .or. size(h) /= size(y)) error stop "fit_sqrt_har_variance: array sizes differ"
        if (ntrain < 30 .or. ntrain > size(y)) error stop "fit_sqrt_har_variance: invalid ntrain"
        if (allocated(obj_y)) deallocate(obj_y)
        if (allocated(obj_x)) deallocate(obj_x)
        allocate(obj_y(size(y)), obj_x(size(x)))
        obj_y = y
        obj_x = max(x, min_var)
        obj_ntrain = ntrain

        mean_y2 = max(sum(y(1:ntrain)**2) / real(ntrain, dp), min_var)
        p = 0.0_dp
        p(1) = sqrt(mean_y2)
        call bfgs_minimize(sqrt_har_obj, p, 4, 300, 1.0e-5_dp, fopt, niter, converged)
        result%coef = p
        result%niter = niter
        result%converged = converged
        call sqrt_har_variance_path(x, result%coef, h)
        result%loglik = gaussian_variance_loglik(y(23:ntrain), h(23:ntrain))
    end subroutine fit_sqrt_har_variance

    subroutine fit_harq_variance(y, x, rq, ntrain, result, h)
        ! Fit HARQ with a demeaned sqrt(realized quarticity) interaction on the daily lag.
        real(dp), intent(in) :: y(:), x(:), rq(:)
        integer, intent(in) :: ntrain
        type(harq_variance_result_t), intent(out) :: result
        real(dp), intent(out) :: h(:)
        real(dp) :: p(5), fopt, mean_y2, mean_x
        integer :: niter
        logical :: converged

        if (size(y) /= size(x) .or. size(y) /= size(rq) .or. size(h) /= size(y)) then
            error stop "fit_harq_variance: array sizes differ"
        end if
        if (ntrain < 30 .or. ntrain > size(y)) error stop "fit_harq_variance: invalid ntrain"
        if (allocated(obj_y)) deallocate(obj_y)
        if (allocated(obj_x)) deallocate(obj_x)
        if (allocated(obj_x2)) deallocate(obj_x2)
        allocate(obj_y(size(y)), obj_x(size(x)), obj_x2(size(rq)))
        obj_y = y
        obj_x = max(x, min_var)
        obj_x2 = max(rq, min_var)
        obj_ntrain = ntrain

        mean_y2 = max(sum(y(1:ntrain)**2) / real(ntrain, dp), min_var)
        mean_x = max(sum(obj_x(1:ntrain)) / real(ntrain, dp), min_var)
        p(1) = log(0.20_dp*mean_y2)
        p(2:4) = log(max(0.20_dp*mean_y2 / mean_x, min_var))
        p(5) = 0.0_dp
        call bfgs_minimize(harq_obj, p, 5, 300, 1.0e-5_dp, fopt, niter, converged)
        result%coef(1:4) = exp(p(1:4))
        result%coef(5) = p(5)
        result%niter = niter
        result%converged = converged
        call harq_variance_path(x, rq, result%coef, ntrain, h)
        result%loglik = gaussian_variance_loglik(y(23:ntrain), h(23:ntrain))
    end subroutine fit_harq_variance

    subroutine fit_harj_variance(y, x, jump, ntrain, result, h)
        ! Fit HARJ with daily, weekly, and monthly realized jump variation terms.
        real(dp), intent(in) :: y(:), x(:), jump(:)
        integer, intent(in) :: ntrain
        type(harj_variance_result_t), intent(out) :: result
        real(dp), intent(out) :: h(:)
        real(dp) :: p(7), fopt, mean_y2, mean_x, mean_jump
        integer :: niter
        logical :: converged

        if (size(y) /= size(x) .or. size(y) /= size(jump) .or. size(h) /= size(y)) then
            error stop "fit_harj_variance: array sizes differ"
        end if
        if (ntrain < 30 .or. ntrain > size(y)) error stop "fit_harj_variance: invalid ntrain"
        if (allocated(obj_y)) deallocate(obj_y)
        if (allocated(obj_x)) deallocate(obj_x)
        if (allocated(obj_x2)) deallocate(obj_x2)
        allocate(obj_y(size(y)), obj_x(size(x)), obj_x2(size(jump)))
        obj_y = y
        obj_x = max(x, min_var)
        obj_x2 = max(jump, 0.0_dp)
        obj_ntrain = ntrain

        mean_y2 = max(sum(y(1:ntrain)**2) / real(ntrain, dp), min_var)
        mean_x = max(sum(obj_x(1:ntrain)) / real(ntrain, dp), min_var)
        mean_jump = max(sum(obj_x2(1:ntrain)) / real(ntrain, dp), min_var)
        p(1) = log(0.15_dp*mean_y2)
        p(2:4) = log(max(0.20_dp*mean_y2 / mean_x, min_var))
        p(5:7) = log(max(0.05_dp*mean_y2 / mean_jump, min_var))
        call bfgs_minimize(harj_obj, p, 7, 300, 1.0e-5_dp, fopt, niter, converged)
        result%coef = exp(p)
        result%niter = niter
        result%converged = converged
        call harj_variance_path(x, jump, result%coef, h)
        result%loglik = gaussian_variance_loglik(y(23:ntrain), h(23:ntrain))
    end subroutine fit_harj_variance

    subroutine fit_harqj_variance(y, x, rq, jump, ntrain, result, h)
        ! Fit HARQJ with HAR, jump, and quarticity interaction terms.
        real(dp), intent(in) :: y(:), x(:), rq(:), jump(:)
        integer, intent(in) :: ntrain
        type(harqj_variance_result_t), intent(out) :: result
        real(dp), intent(out) :: h(:)
        real(dp) :: p(8), fopt, mean_y2, mean_x, mean_jump
        integer :: niter
        logical :: converged

        if (size(y) /= size(x) .or. size(y) /= size(rq) .or. size(y) /= size(jump) .or. size(h) /= size(y)) then
            error stop "fit_harqj_variance: array sizes differ"
        end if
        if (ntrain < 30 .or. ntrain > size(y)) error stop "fit_harqj_variance: invalid ntrain"
        if (allocated(obj_y)) deallocate(obj_y)
        if (allocated(obj_x)) deallocate(obj_x)
        if (allocated(obj_x2)) deallocate(obj_x2)
        if (allocated(obj_x3)) deallocate(obj_x3)
        allocate(obj_y(size(y)), obj_x(size(x)), obj_x2(size(rq)), obj_x3(size(jump)))
        obj_y = y
        obj_x = max(x, min_var)
        obj_x2 = max(rq, min_var)
        obj_x3 = max(jump, 0.0_dp)
        obj_ntrain = ntrain

        mean_y2 = max(sum(y(1:ntrain)**2) / real(ntrain, dp), min_var)
        mean_x = max(sum(obj_x(1:ntrain)) / real(ntrain, dp), min_var)
        mean_jump = max(sum(obj_x3(1:ntrain)) / real(ntrain, dp), min_var)
        p(1) = log(0.15_dp*mean_y2)
        p(2:4) = log(max(0.20_dp*mean_y2 / mean_x, min_var))
        p(5:7) = log(max(0.05_dp*mean_y2 / mean_jump, min_var))
        p(8) = 0.0_dp
        call bfgs_minimize(harqj_obj, p, 8, 300, 1.0e-5_dp, fopt, niter, converged)
        result%coef(1:7) = exp(p(1:7))
        result%coef(8) = p(8)
        result%niter = niter
        result%converged = converged
        call harqj_variance_path(x, rq, jump, result%coef, ntrain, h)
        result%loglik = gaussian_variance_loglik(y(23:ntrain), h(23:ntrain))
    end subroutine fit_harqj_variance

    subroutine fit_har_negret_variance(y, x, ntrain, result, h)
        ! Fit positive HAR plus lagged negative-return-squared term.
        real(dp), intent(in) :: y(:), x(:)
        integer, intent(in) :: ntrain
        type(har_negret_variance_result_t), intent(out) :: result
        real(dp), intent(out) :: h(:)
        real(dp) :: p(5), fopt, mean_y2, mean_x
        integer :: niter
        logical :: converged

        if (size(y) /= size(x) .or. size(h) /= size(y)) error stop "fit_har_negret_variance: array sizes differ"
        if (ntrain < 30 .or. ntrain > size(y)) error stop "fit_har_negret_variance: invalid ntrain"
        if (allocated(obj_y)) deallocate(obj_y)
        if (allocated(obj_x)) deallocate(obj_x)
        allocate(obj_y(size(y)), obj_x(size(x)))
        obj_y = y
        obj_x = max(x, min_var)
        obj_ntrain = ntrain

        mean_y2 = max(sum(y(1:ntrain)**2) / real(ntrain, dp), min_var)
        mean_x = max(sum(obj_x(1:ntrain)) / real(ntrain, dp), min_var)
        p(1) = log(0.20_dp*mean_y2)
        p(2:4) = log(max(0.20_dp*mean_y2 / mean_x, min_var))
        p(5) = log(0.20_dp)
        call bfgs_minimize(har_negret_obj, p, 5, 300, 1.0e-5_dp, fopt, niter, converged)
        result%coef = exp(p)
        result%niter = niter
        result%converged = converged
        call har_negret_variance_path(y, x, result%coef, h)
        result%loglik = gaussian_variance_loglik(y(23:ntrain), h(23:ntrain))
    end subroutine fit_har_negret_variance

    subroutine fit_har_lev_variance(y, x, ntrain, result, h)
        ! Fit HAR plus daily, weekly, and monthly negative-return-squared leverage terms.
        real(dp), intent(in) :: y(:), x(:)
        integer, intent(in) :: ntrain
        type(har_lev_variance_result_t), intent(out) :: result
        real(dp), intent(out) :: h(:)
        real(dp) :: p(7), fopt, mean_y2, mean_x, mean_neg
        integer :: niter
        logical :: converged

        if (size(y) /= size(x) .or. size(h) /= size(y)) error stop "fit_har_lev_variance: array sizes differ"
        if (ntrain < 30 .or. ntrain > size(y)) error stop "fit_har_lev_variance: invalid ntrain"
        if (allocated(obj_y)) deallocate(obj_y)
        if (allocated(obj_x)) deallocate(obj_x)
        allocate(obj_y(size(y)), obj_x(size(x)))
        obj_y = y
        obj_x = max(x, min_var)
        obj_ntrain = ntrain

        mean_y2 = max(sum(y(1:ntrain)**2) / real(ntrain, dp), min_var)
        mean_x = max(sum(obj_x(1:ntrain)) / real(ntrain, dp), min_var)
        mean_neg = max(sum(min(y(1:ntrain), 0.0_dp)**2) / real(ntrain, dp), min_var)
        p(1) = log(0.15_dp*mean_y2)
        p(2:4) = log(max(0.15_dp*mean_y2 / mean_x, min_var))
        p(5:7) = log(max(0.15_dp*mean_y2 / mean_neg, min_var))
        call bfgs_minimize(har_lev_obj, p, 7, 300, 1.0e-5_dp, fopt, niter, converged)
        result%coef = exp(p)
        result%niter = niter
        result%converged = converged
        call har_lev_variance_path(y, x, result%coef, h)
        result%loglik = gaussian_variance_loglik(y(23:ntrain), h(23:ntrain))
    end subroutine fit_har_lev_variance

    subroutine fit_semivar_har_variance(y, x_pos, x_neg, ntrain, result, h)
        ! Fit positive-coefficient HAR using separate positive and negative semivariance predictors.
        real(dp), intent(in) :: y(:), x_pos(:), x_neg(:)
        integer, intent(in) :: ntrain
        type(semivar_har_variance_result_t), intent(out) :: result
        real(dp), intent(out) :: h(:)
        real(dp) :: p(7), fopt, mean_y2, mean_x
        integer :: niter
        logical :: converged

        if (size(y) /= size(x_pos) .or. size(y) /= size(x_neg) .or. size(y) /= size(h)) then
            error stop "fit_semivar_har_variance: array sizes differ"
        end if
        if (ntrain < 30 .or. ntrain > size(y)) error stop "fit_semivar_har_variance: invalid ntrain"
        if (allocated(obj_y)) deallocate(obj_y)
        if (allocated(obj_x)) deallocate(obj_x)
        if (allocated(obj_x2)) deallocate(obj_x2)
        allocate(obj_y(size(y)), obj_x(size(x_pos)), obj_x2(size(x_neg)))
        obj_y = y
        obj_x = max(x_pos, min_var)
        obj_x2 = max(x_neg, min_var)
        obj_ntrain = ntrain

        mean_y2 = max(sum(y(1:ntrain)**2) / real(ntrain, dp), min_var)
        mean_x = max(sum(obj_x(1:ntrain) + obj_x2(1:ntrain)) / real(ntrain, dp), min_var)
        p(1) = log(0.15_dp*mean_y2)
        p(2:7) = log(max(0.15_dp*mean_y2 / mean_x, min_var))
        call bfgs_minimize(semivar_har_obj, p, 7, 300, 1.0e-5_dp, fopt, niter, converged)
        result%coef = exp(p)
        result%niter = niter
        result%converged = converged
        call semivar_har_variance_path(x_pos, x_neg, result%coef, h)
        result%loglik = gaussian_variance_loglik(y(23:ntrain), h(23:ntrain))
    end subroutine fit_semivar_har_variance

    subroutine fit_midas_variance(y, x, ntrain, k_lag, result, h)
        ! Fit h_t = a + b*sum_k w_k(theta1,theta2)*x_{t-k} by Gaussian likelihood.
        real(dp), intent(in) :: y(:), x(:)
        integer, intent(in) :: ntrain, k_lag
        type(midas_variance_result_t), intent(out) :: result
        real(dp), intent(out) :: h(:)
        real(dp) :: p(4), fopt, mean_y2, mean_x
        integer :: niter
        logical :: converged

        if (size(y) /= size(x) .or. size(y) /= size(h)) error stop "fit_midas_variance: array sizes differ"
        if (k_lag < 2) error stop "fit_midas_variance: k_lag must be at least 2"
        if (ntrain <= k_lag + 5 .or. ntrain > size(y)) error stop "fit_midas_variance: invalid ntrain"
        if (allocated(obj_y)) deallocate(obj_y)
        if (allocated(obj_x)) deallocate(obj_x)
        allocate(obj_y(size(y)), obj_x(size(x)))
        obj_y = y
        obj_x = max(x, min_var)
        obj_ntrain = ntrain
        obj_k_lag = k_lag

        mean_y2 = max(sum(y(1:ntrain)**2) / real(ntrain, dp), min_var)
        mean_x = max(sum(obj_x(1:ntrain)) / real(ntrain, dp), min_var)
        p(1) = log(0.1_dp*mean_y2)
        p(2) = log(max(0.9_dp*mean_y2 / mean_x, min_var))
        p(3) = log(1.0_dp)
        p(4) = log(2.0_dp)
        call bfgs_minimize(midas_obj, p, 4, 300, 1.0e-5_dp, fopt, niter, converged)
        result%a = exp(p(1))
        result%b = exp(p(2))
        result%theta1 = exp(p(3))
        result%theta2 = exp(p(4))
        result%niter = niter
        result%converged = converged
        call midas_variance_path(x, k_lag, result%a, result%b, result%theta1, result%theta2, h)
        result%loglik = gaussian_variance_loglik(y(k_lag + 1:ntrain), h(k_lag + 1:ntrain))
    end subroutine fit_midas_variance

    subroutine fit_heavy_variance(y, x, ntrain, max_iter, gtol, result, h)
        ! Fit the HEAVY return and realized-measure variance equations.
        real(dp), intent(in) :: y(:), x(:), gtol
        integer, intent(in) :: ntrain, max_iter
        type(heavy_variance_result_t), intent(out) :: result
        real(dp), intent(out) :: h(:)
        real(dp) :: p_ret(3), p_rm(3), grad(3), fopt_ret, fopt_rm
        real(dp) :: mu_y2, mu_x
        integer :: niter_ret, niter_rm
        logical :: conv_ret, conv_rm

        if (size(y) /= size(x) .or. size(h) /= size(y)) error stop "fit_heavy_variance: array sizes differ"
        if (ntrain < 5 .or. ntrain > size(y)) error stop "fit_heavy_variance: invalid ntrain"
        if (max_iter < 1) error stop "fit_heavy_variance: max_iter must be positive"
        if (allocated(obj_y)) deallocate(obj_y)
        if (allocated(obj_x)) deallocate(obj_x)
        allocate(obj_y(size(y)), obj_x(size(x)))
        obj_y = y
        obj_x = max(x, min_var)
        obj_ntrain = ntrain

        mu_y2 = max(sum(y(1:ntrain)**2) / real(ntrain, dp), min_var)
        mu_x = max(sum(obj_x(1:ntrain)) / real(ntrain, dp), min_var)
        call heavy_inverse_transform(max(0.1_dp*mu_y2, min_var), 0.30_dp*mu_y2/mu_x, 0.50_dp, p_ret)
        call heavy_inverse_transform(max(0.1_dp*mu_x, min_var), 0.60_dp, 0.30_dp, p_rm)

        call bfgs_minimize(heavy_return_obj, p_ret, 3, max_iter, gtol, fopt_ret, niter_ret, conv_ret)
        call bfgs_minimize(heavy_rm_obj, p_rm, 3, max_iter, gtol, fopt_rm, niter_rm, conv_rm)
        call heavy_transform(p_ret, result%omega, result%alpha, result%beta)
        call heavy_transform(p_rm, result%omega_rm, result%alpha_rm, result%beta_rm)
        call heavy_variance_path(x, result%omega, result%alpha, result%beta, ntrain, h)
        result%loglik = gaussian_variance_loglik(y(2:ntrain), h(2:ntrain))
        result%loglik_rm = heavy_rm_loglik(result%omega_rm, result%alpha_rm, result%beta_rm)
        result%niter = niter_ret + niter_rm
        result%converged = conv_ret .and. conv_rm
        call heavy_return_obj(p_ret, 3, fopt_ret, grad)
        call heavy_rm_obj(p_rm, 3, fopt_rm, grad)
    end subroutine fit_heavy_variance

    subroutine har_variance_path(x, coef, h)
        ! Compute HAR variance forecasts from a realized variance measure.
        real(dp), intent(in) :: x(:), coef(4)
        real(dp), intent(out) :: h(:)
        real(dp) :: row(4), fallback
        integer :: t

        if (size(x) /= size(h)) error stop "har_variance_path: array sizes differ"
        fallback = max(sum(max(x, min_var)) / real(size(x), dp), min_var)
        h(1:min(22, size(h))) = fallback
        do t = 23, size(x)
            call har_row(x, t, row)
            h(t) = max(sum(coef*row), min_var)
        end do
    end subroutine har_variance_path

    subroutine harx_variance_path(x, x_exog, coef, h)
        ! Compute HAR forecasts plus an exogenous variance predictor known at forecast time.
        real(dp), intent(in) :: x(:), x_exog(:), coef(5)
        real(dp), intent(out) :: h(:)
        real(dp) :: row(5), fallback
        integer :: t

        if (size(x) /= size(x_exog) .or. size(x) /= size(h)) error stop "harx_variance_path: array sizes differ"
        fallback = max(sum(max(x, min_var)) / real(size(x), dp), min_var)
        h(1:min(22, size(h))) = fallback
        do t = 23, size(x)
            call harx_row(x, x_exog, t, row)
            h(t) = max(sum(coef*row), min_var)
        end do
    end subroutine harx_variance_path

    subroutine harx_lev_variance_path(y, x, x_exog, coef, h)
        ! Compute HARX forecasts plus daily/weekly/monthly negative-return-squared terms.
        real(dp), intent(in) :: y(:), x(:), x_exog(:), coef(8)
        real(dp), intent(out) :: h(:)
        real(dp) :: row(8), fallback
        integer :: t

        if (size(y) /= size(x) .or. size(y) /= size(x_exog) .or. size(y) /= size(h)) then
            error stop "harx_lev_variance_path: array sizes differ"
        end if
        fallback = max(sum(max(x, min_var)) / real(size(x), dp), min_var)
        h(1:min(22, size(h))) = fallback
        do t = 23, size(x)
            call harx_lev_row(y, x, x_exog, t, row)
            h(t) = max(sum(coef*row), min_var)
        end do
    end subroutine harx_lev_variance_path

    subroutine log_har_variance_path(x, coef, h)
        ! Compute log-HAR variance forecasts from a realized variance measure.
        real(dp), intent(in) :: x(:), coef(4)
        real(dp), intent(out) :: h(:)
        real(dp) :: row(4), fallback, logh
        integer :: t

        if (size(x) /= size(h)) error stop "log_har_variance_path: array sizes differ"
        fallback = max(sum(max(x, min_var)) / real(size(x), dp), min_var)
        h(1:min(22, size(h))) = fallback
        do t = 23, size(x)
            call log_har_row(x, t, row)
            logh = sum(coef*row)
            h(t) = max(exp(min(max(logh, log(min_var)), 50.0_dp)), min_var)
        end do
    end subroutine log_har_variance_path

    subroutine sqrt_har_variance_path(x, coef, h)
        ! Compute sqrt-HAR variance forecasts from a realized variance measure.
        real(dp), intent(in) :: x(:), coef(4)
        real(dp), intent(out) :: h(:)
        real(dp) :: row(4), fallback, vol
        integer :: t

        if (size(x) /= size(h)) error stop "sqrt_har_variance_path: array sizes differ"
        fallback = max(sum(max(x, min_var)) / real(size(x), dp), min_var)
        h(1:min(22, size(h))) = fallback
        do t = 23, size(x)
            call sqrt_har_row(x, t, row)
            vol = sum(coef*row)
            h(t) = max(vol**2, min_var)
        end do
    end subroutine sqrt_har_variance_path

    subroutine harq_variance_path(x, rq, coef, ninit, h)
        ! Compute HARQ forecasts with a one-day quarticity interaction term.
        real(dp), intent(in) :: x(:), rq(:), coef(5)
        integer, intent(in) :: ninit
        real(dp), intent(out) :: h(:)
        real(dp) :: row(4), fallback, rq_center, rq_term
        integer :: t

        if (size(x) /= size(rq) .or. size(x) /= size(h)) error stop "harq_variance_path: array sizes differ"
        if (ninit < 23 .or. ninit > size(x)) error stop "harq_variance_path: invalid ninit"
        fallback = max(sum(max(x, min_var)) / real(size(x), dp), min_var)
        rq_center = sqrt(max(sum(max(rq(1:ninit), min_var)) / real(ninit, dp), min_var))
        h(1:min(22, size(h))) = fallback
        do t = 23, size(x)
            call har_row(x, t, row)
            rq_term = (sqrt(max(rq(t - 1), min_var)) - rq_center) * row(2)
            h(t) = max(sum(coef(1:4)*row) + coef(5)*rq_term, min_var)
        end do
    end subroutine harq_variance_path

    subroutine harj_variance_path(x, jump, coef, h)
        ! Compute HARJ forecasts from HAR variance terms plus realized jump variation.
        real(dp), intent(in) :: x(:), jump(:), coef(7)
        real(dp), intent(out) :: h(:)
        real(dp) :: row(7), fallback
        integer :: t

        if (size(x) /= size(jump) .or. size(x) /= size(h)) error stop "harj_variance_path: array sizes differ"
        fallback = max(sum(max(x, min_var)) / real(size(x), dp), min_var)
        h(1:min(22, size(h))) = fallback
        do t = 23, size(x)
            call harj_row(x, jump, t, row)
            h(t) = max(sum(coef*row), min_var)
        end do
    end subroutine harj_variance_path

    subroutine harqj_variance_path(x, rq, jump, coef, ninit, h)
        ! Compute HARQJ forecasts from HAR, jump, and quarticity interaction terms.
        real(dp), intent(in) :: x(:), rq(:), jump(:), coef(8)
        integer, intent(in) :: ninit
        real(dp), intent(out) :: h(:)
        real(dp) :: row(7), fallback, rq_center, rq_term
        integer :: t

        if (size(x) /= size(rq) .or. size(x) /= size(jump) .or. size(x) /= size(h)) then
            error stop "harqj_variance_path: array sizes differ"
        end if
        if (ninit < 23 .or. ninit > size(x)) error stop "harqj_variance_path: invalid ninit"
        fallback = max(sum(max(x, min_var)) / real(size(x), dp), min_var)
        rq_center = sqrt(max(sum(max(rq(1:ninit), min_var)) / real(ninit, dp), min_var))
        h(1:min(22, size(h))) = fallback
        do t = 23, size(x)
            call harj_row(x, jump, t, row)
            rq_term = (sqrt(max(rq(t - 1), min_var)) - rq_center) * row(2)
            h(t) = max(sum(coef(1:7)*row) + coef(8)*rq_term, min_var)
        end do
    end subroutine harqj_variance_path

    subroutine har_negret_variance_path(y, x, coef, h)
        ! Compute HAR variance forecasts plus lagged negative daily return squared.
        real(dp), intent(in) :: y(:), x(:), coef(5)
        real(dp), intent(out) :: h(:)
        real(dp) :: row(4), fallback, negsq
        integer :: t

        if (size(y) /= size(x) .or. size(y) /= size(h)) error stop "har_negret_variance_path: array sizes differ"
        fallback = max(sum(max(x, min_var)) / real(size(x), dp), min_var)
        h(1:min(22, size(h))) = fallback
        do t = 23, size(x)
            call har_row(x, t, row)
            negsq = min(y(t - 1), 0.0_dp)**2
            h(t) = max(sum(coef(1:4)*row) + coef(5)*negsq, min_var)
        end do
    end subroutine har_negret_variance_path

    subroutine har_lev_variance_path(y, x, coef, h)
        ! Compute HAR forecasts plus daily/weekly/monthly negative-return-squared terms.
        real(dp), intent(in) :: y(:), x(:), coef(7)
        real(dp), intent(out) :: h(:)
        real(dp) :: row(7), fallback
        integer :: t

        if (size(y) /= size(x) .or. size(y) /= size(h)) error stop "har_lev_variance_path: array sizes differ"
        fallback = max(sum(max(x, min_var)) / real(size(x), dp), min_var)
        h(1:min(22, size(h))) = fallback
        do t = 23, size(x)
            call har_lev_row(y, x, t, row)
            h(t) = max(sum(coef*row), min_var)
        end do
    end subroutine har_lev_variance_path

    subroutine semivar_har_variance_path(x_pos, x_neg, coef, h)
        ! Compute semivariance-HAR variance forecasts.
        real(dp), intent(in) :: x_pos(:), x_neg(:), coef(7)
        real(dp), intent(out) :: h(:)
        real(dp) :: row(7), fallback
        integer :: t

        if (size(x_pos) /= size(x_neg) .or. size(x_pos) /= size(h)) then
            error stop "semivar_har_variance_path: array sizes differ"
        end if
        fallback = max(sum(max(x_pos + x_neg, min_var)) / real(size(x_pos), dp), min_var)
        h(1:min(22, size(h))) = fallback
        do t = 23, size(x_pos)
            call semivar_har_row(x_pos, x_neg, t, row)
            h(t) = max(sum(coef*row), min_var)
        end do
    end subroutine semivar_har_variance_path

    subroutine midas_variance_path(x, k_lag, a, bcoef, theta1, theta2, h)
        ! Compute positive Beta-weight MIDAS variance forecasts.
        real(dp), intent(in) :: x(:), a, bcoef, theta1, theta2
        integer, intent(in) :: k_lag
        real(dp), intent(out) :: h(:)
        real(dp), allocatable :: weights(:)
        real(dp) :: xw, fallback
        integer :: t, k

        if (size(x) /= size(h)) error stop "midas_variance_path: array sizes differ"
        if (k_lag < 2 .or. k_lag >= size(x)) error stop "midas_variance_path: invalid k_lag"
        allocate(weights(k_lag))
        call beta_lag_weights(k_lag, theta1, theta2, weights)
        fallback = max(sum(max(x, min_var)) / real(size(x), dp), min_var)
        h(1:k_lag) = max(a + bcoef*fallback, min_var)
        do t = k_lag + 1, size(x)
            xw = 0.0_dp
            do k = 1, k_lag
                xw = xw + weights(k)*max(x(t - k), min_var)
            end do
            h(t) = max(a + bcoef*xw, min_var)
        end do
        deallocate(weights)
    end subroutine midas_variance_path

    subroutine heavy_variance_path(x, omega, alpha, beta, ninit, h)
        ! Compute HEAVY h_t = omega + alpha*x_{t-1} + beta*h_{t-1}.
        real(dp), intent(in) :: x(:), omega, alpha, beta
        integer, intent(in) :: ninit
        real(dp), intent(out) :: h(:)
        integer :: t

        if (size(x) /= size(h)) error stop "heavy_variance_path: array sizes differ"
        if (ninit < 1 .or. ninit > size(x)) error stop "heavy_variance_path: invalid ninit"
        h(1) = max(sum(max(x(1:ninit), min_var)) / real(ninit, dp), min_var)
        do t = 2, size(x)
            h(t) = max(omega + alpha*max(x(t - 1), min_var) + beta*h(t - 1), min_var)
        end do
    end subroutine heavy_variance_path

    subroutine affine_variance_path(x, a, b, h)
        ! Compute h_t = a + b*x_t.
        real(dp), intent(in) :: x(:), a, b
        real(dp), intent(out) :: h(:)

        if (size(x) /= size(h)) error stop "affine_variance_path: array sizes differ"
        h = max(a + b*max(x, min_var), min_var)
    end subroutine affine_variance_path

    subroutine affine2_variance_path(x1, x2, a, b1, b2, h)
        ! Compute h_t = a + b1*x1_t + b2*x2_t.
        real(dp), intent(in) :: x1(:), x2(:), a, b1, b2
        real(dp), intent(out) :: h(:)

        if (size(x1) /= size(x2) .or. size(x1) /= size(h)) error stop "affine2_variance_path: array sizes differ"
        h = max(a + b1*max(x1, min_var) + b2*max(x2, min_var), min_var)
    end subroutine affine2_variance_path

    subroutine ewma_affine_variance_path(x, lambda, a, b, ninit, h)
        ! Compute h_t = a + b*EWMA(x_{t-1}) using mean x(1:ninit) as startup.
        real(dp), intent(in) :: x(:), lambda, a, b
        integer, intent(in) :: ninit
        real(dp), intent(out) :: h(:)
        real(dp) :: lam, m
        integer :: t

        if (size(x) /= size(h)) error stop "ewma_affine_variance_path: array sizes differ"
        if (ninit < 1 .or. ninit > size(x)) error stop "ewma_affine_variance_path: invalid ninit"
        lam = min(max(lambda, 0.0001_dp), 0.9999_dp)
        m = max(sum(max(x(1:ninit), min_var)) / real(ninit, dp), min_var)
        h(1) = max(a + b*m, min_var)
        do t = 2, size(x)
            m = max(lam*m + (1.0_dp - lam)*max(x(t - 1), min_var), min_var)
            h(t) = max(a + b*m, min_var)
        end do
    end subroutine ewma_affine_variance_path

    real(dp) function gaussian_variance_loglik(y, h)
        ! Gaussian zero-mean log likelihood for return y and variance h.
        real(dp), intent(in) :: y(:), h(:)
        integer :: i

        if (size(y) /= size(h)) error stop "gaussian_variance_loglik: array sizes differ"
        gaussian_variance_loglik = 0.0_dp
        do i = 1, size(y)
            gaussian_variance_loglik = gaussian_variance_loglik - log_sqrt_2pi - 0.5_dp*log(max(h(i), min_var)) - &
                                       0.5_dp*y(i)**2 / max(h(i), min_var)
        end do
    end function gaussian_variance_loglik

    real(dp) function qlike_loss(y, h)
        ! Average QLIKE loss log(h_t) + y_t^2/h_t.
        real(dp), intent(in) :: y(:), h(:)

        if (size(y) /= size(h)) error stop "qlike_loss: array sizes differ"
        qlike_loss = sum(log(max(h, min_var)) + y**2 / max(h, min_var)) / real(size(y), dp)
    end function qlike_loss

    subroutine ewma_affine_obj(p, np, f, g)
        integer, intent(in) :: np
        real(dp), intent(in) :: p(np)
        real(dp), intent(out) :: f, g(np)
        real(dp), allocatable :: pp(:), pm(:)
        real(dp) :: step, fp, fm
        integer :: j

        f = ewma_affine_nll(p)
        allocate(pp(np), pm(np))
        do j = 1, np
            step = 1.0e-5_dp * max(1.0_dp, abs(p(j)))
            pp = p
            pm = p
            pp(j) = pp(j) + step
            pm(j) = pm(j) - step
            fp = ewma_affine_nll(pp)
            fm = ewma_affine_nll(pm)
            g(j) = (fp - fm) / (2.0_dp*step)
        end do
        deallocate(pp, pm)
    end subroutine ewma_affine_obj

    subroutine affine_obj(p, np, f, g)
        integer, intent(in) :: np
        real(dp), intent(in) :: p(np)
        real(dp), intent(out) :: f, g(np)

        call finite_diff_obj(affine_nll, p, np, f, g)
    end subroutine affine_obj

    subroutine affine2_obj(p, np, f, g)
        integer, intent(in) :: np
        real(dp), intent(in) :: p(np)
        real(dp), intent(out) :: f, g(np)

        call finite_diff_obj(affine2_nll, p, np, f, g)
    end subroutine affine2_obj

    subroutine har_obj(p, np, f, g)
        integer, intent(in) :: np
        real(dp), intent(in) :: p(np)
        real(dp), intent(out) :: f, g(np)
        real(dp), allocatable :: pp(:), pm(:)
        real(dp) :: step, fp, fm
        integer :: j

        f = har_nll(p)
        allocate(pp(np), pm(np))
        do j = 1, np
            step = 1.0e-5_dp * max(1.0_dp, abs(p(j)))
            pp = p
            pm = p
            pp(j) = pp(j) + step
            pm(j) = pm(j) - step
            fp = har_nll(pp)
            fm = har_nll(pm)
            g(j) = (fp - fm) / (2.0_dp*step)
        end do
        deallocate(pp, pm)
    end subroutine har_obj

    subroutine harx_obj(p, np, f, g)
        integer, intent(in) :: np
        real(dp), intent(in) :: p(np)
        real(dp), intent(out) :: f, g(np)

        call finite_diff_obj(harx_nll, p, np, f, g)
    end subroutine harx_obj

    subroutine harx_lev_obj(p, np, f, g)
        integer, intent(in) :: np
        real(dp), intent(in) :: p(np)
        real(dp), intent(out) :: f, g(np)

        call finite_diff_obj(harx_lev_nll, p, np, f, g)
    end subroutine harx_lev_obj

    subroutine log_har_obj(p, np, f, g)
        integer, intent(in) :: np
        real(dp), intent(in) :: p(np)
        real(dp), intent(out) :: f, g(np)

        call finite_diff_obj(log_har_nll, p, np, f, g)
    end subroutine log_har_obj

    subroutine sqrt_har_obj(p, np, f, g)
        integer, intent(in) :: np
        real(dp), intent(in) :: p(np)
        real(dp), intent(out) :: f, g(np)

        call finite_diff_obj(sqrt_har_nll, p, np, f, g)
    end subroutine sqrt_har_obj

    subroutine harq_obj(p, np, f, g)
        integer, intent(in) :: np
        real(dp), intent(in) :: p(np)
        real(dp), intent(out) :: f, g(np)

        call finite_diff_obj(harq_nll, p, np, f, g)
    end subroutine harq_obj

    subroutine harj_obj(p, np, f, g)
        integer, intent(in) :: np
        real(dp), intent(in) :: p(np)
        real(dp), intent(out) :: f, g(np)

        call finite_diff_obj(harj_nll, p, np, f, g)
    end subroutine harj_obj

    subroutine harqj_obj(p, np, f, g)
        integer, intent(in) :: np
        real(dp), intent(in) :: p(np)
        real(dp), intent(out) :: f, g(np)

        call finite_diff_obj(harqj_nll, p, np, f, g)
    end subroutine harqj_obj

    subroutine har_negret_obj(p, np, f, g)
        integer, intent(in) :: np
        real(dp), intent(in) :: p(np)
        real(dp), intent(out) :: f, g(np)
        real(dp), allocatable :: pp(:), pm(:)
        real(dp) :: step, fp, fm
        integer :: j

        f = har_negret_nll(p)
        allocate(pp(np), pm(np))
        do j = 1, np
            step = 1.0e-5_dp * max(1.0_dp, abs(p(j)))
            pp = p
            pm = p
            pp(j) = pp(j) + step
            pm(j) = pm(j) - step
            fp = har_negret_nll(pp)
            fm = har_negret_nll(pm)
            g(j) = (fp - fm) / (2.0_dp*step)
        end do
        deallocate(pp, pm)
    end subroutine har_negret_obj

    subroutine har_lev_obj(p, np, f, g)
        integer, intent(in) :: np
        real(dp), intent(in) :: p(np)
        real(dp), intent(out) :: f, g(np)

        call finite_diff_obj(har_lev_nll, p, np, f, g)
    end subroutine har_lev_obj

    subroutine semivar_har_obj(p, np, f, g)
        integer, intent(in) :: np
        real(dp), intent(in) :: p(np)
        real(dp), intent(out) :: f, g(np)
        real(dp), allocatable :: pp(:), pm(:)
        real(dp) :: step, fp, fm
        integer :: j

        f = semivar_har_nll(p)
        allocate(pp(np), pm(np))
        do j = 1, np
            step = 1.0e-5_dp * max(1.0_dp, abs(p(j)))
            pp = p
            pm = p
            pp(j) = pp(j) + step
            pm(j) = pm(j) - step
            fp = semivar_har_nll(pp)
            fm = semivar_har_nll(pm)
            g(j) = (fp - fm) / (2.0_dp*step)
        end do
        deallocate(pp, pm)
    end subroutine semivar_har_obj

    subroutine midas_obj(p, np, f, g)
        integer, intent(in) :: np
        real(dp), intent(in) :: p(np)
        real(dp), intent(out) :: f, g(np)
        real(dp), allocatable :: pp(:), pm(:)
        real(dp) :: step, fp, fm
        integer :: j

        f = midas_nll(p)
        allocate(pp(np), pm(np))
        do j = 1, np
            step = 1.0e-5_dp * max(1.0_dp, abs(p(j)))
            pp = p
            pm = p
            pp(j) = pp(j) + step
            pm(j) = pm(j) - step
            fp = midas_nll(pp)
            fm = midas_nll(pm)
            g(j) = (fp - fm) / (2.0_dp*step)
        end do
        deallocate(pp, pm)
    end subroutine midas_obj

    subroutine heavy_return_obj(p, np, f, g)
        integer, intent(in) :: np
        real(dp), intent(in) :: p(np)
        real(dp), intent(out) :: f, g(np)

        call finite_diff_obj(heavy_return_nll, p, np, f, g)
    end subroutine heavy_return_obj

    subroutine heavy_rm_obj(p, np, f, g)
        integer, intent(in) :: np
        real(dp), intent(in) :: p(np)
        real(dp), intent(out) :: f, g(np)

        call finite_diff_obj(heavy_rm_nll, p, np, f, g)
    end subroutine heavy_rm_obj

    real(dp) function har_nll(p)
        real(dp), intent(in) :: p(:)
        real(dp), allocatable :: h(:)
        real(dp) :: coef(4)

        coef = exp(p)
        allocate(h(size(obj_y)))
        call har_variance_path(obj_x, coef, h)
        har_nll = -gaussian_variance_loglik(obj_y(23:obj_ntrain), h(23:obj_ntrain)) / real(obj_ntrain - 22, dp)
        if (har_nll /= har_nll .or. har_nll > 1.0e29_dp) har_nll = 1.0e30_dp
        deallocate(h)
    end function har_nll

    real(dp) function harx_nll(p)
        real(dp), intent(in) :: p(:)
        real(dp), allocatable :: h(:)
        real(dp) :: coef(5)

        coef = exp(p)
        allocate(h(size(obj_y)))
        call harx_variance_path(obj_x, obj_x2, coef, h)
        harx_nll = -gaussian_variance_loglik(obj_y(23:obj_ntrain), h(23:obj_ntrain)) / real(obj_ntrain - 22, dp)
        if (harx_nll /= harx_nll .or. harx_nll > 1.0e29_dp) harx_nll = 1.0e30_dp
        deallocate(h)
    end function harx_nll

    real(dp) function harx_lev_nll(p)
        real(dp), intent(in) :: p(:)
        real(dp), allocatable :: h(:)
        real(dp) :: coef(8)

        coef = exp(p)
        allocate(h(size(obj_y)))
        call harx_lev_variance_path(obj_y, obj_x, obj_x2, coef, h)
        harx_lev_nll = -gaussian_variance_loglik(obj_y(23:obj_ntrain), h(23:obj_ntrain)) / &
                       real(obj_ntrain - 22, dp)
        if (harx_lev_nll /= harx_lev_nll .or. harx_lev_nll > 1.0e29_dp) harx_lev_nll = 1.0e30_dp
        deallocate(h)
    end function harx_lev_nll

    real(dp) function log_har_nll(p)
        real(dp), intent(in) :: p(:)
        real(dp), allocatable :: h(:)

        allocate(h(size(obj_y)))
        call log_har_variance_path(obj_x, p(1:4), h)
        log_har_nll = -gaussian_variance_loglik(obj_y(23:obj_ntrain), h(23:obj_ntrain)) / real(obj_ntrain - 22, dp)
        if (log_har_nll /= log_har_nll .or. log_har_nll > 1.0e29_dp) log_har_nll = 1.0e30_dp
        deallocate(h)
    end function log_har_nll

    real(dp) function sqrt_har_nll(p)
        real(dp), intent(in) :: p(:)
        real(dp), allocatable :: h(:)

        allocate(h(size(obj_y)))
        call sqrt_har_variance_path(obj_x, p(1:4), h)
        sqrt_har_nll = -gaussian_variance_loglik(obj_y(23:obj_ntrain), h(23:obj_ntrain)) / real(obj_ntrain - 22, dp)
        if (sqrt_har_nll /= sqrt_har_nll .or. sqrt_har_nll > 1.0e29_dp) sqrt_har_nll = 1.0e30_dp
        deallocate(h)
    end function sqrt_har_nll

    real(dp) function harq_nll(p)
        real(dp), intent(in) :: p(:)
        real(dp), allocatable :: h(:)
        real(dp) :: coef(5)

        coef(1:4) = exp(p(1:4))
        coef(5) = p(5)
        allocate(h(size(obj_y)))
        call harq_variance_path(obj_x, obj_x2, coef, obj_ntrain, h)
        harq_nll = -gaussian_variance_loglik(obj_y(23:obj_ntrain), h(23:obj_ntrain)) / real(obj_ntrain - 22, dp)
        if (harq_nll /= harq_nll .or. harq_nll > 1.0e29_dp) harq_nll = 1.0e30_dp
        deallocate(h)
    end function harq_nll

    real(dp) function harj_nll(p)
        real(dp), intent(in) :: p(:)
        real(dp), allocatable :: h(:)
        real(dp) :: coef(7)

        coef = exp(p)
        allocate(h(size(obj_y)))
        call harj_variance_path(obj_x, obj_x2, coef, h)
        harj_nll = -gaussian_variance_loglik(obj_y(23:obj_ntrain), h(23:obj_ntrain)) / real(obj_ntrain - 22, dp)
        if (harj_nll /= harj_nll .or. harj_nll > 1.0e29_dp) harj_nll = 1.0e30_dp
        deallocate(h)
    end function harj_nll

    real(dp) function harqj_nll(p)
        real(dp), intent(in) :: p(:)
        real(dp), allocatable :: h(:)
        real(dp) :: coef(8)

        coef(1:7) = exp(p(1:7))
        coef(8) = p(8)
        allocate(h(size(obj_y)))
        call harqj_variance_path(obj_x, obj_x2, obj_x3, coef, obj_ntrain, h)
        harqj_nll = -gaussian_variance_loglik(obj_y(23:obj_ntrain), h(23:obj_ntrain)) / real(obj_ntrain - 22, dp)
        if (harqj_nll /= harqj_nll .or. harqj_nll > 1.0e29_dp) harqj_nll = 1.0e30_dp
        deallocate(h)
    end function harqj_nll

    real(dp) function har_negret_nll(p)
        real(dp), intent(in) :: p(:)
        real(dp), allocatable :: h(:)
        real(dp) :: coef(5)

        coef = exp(p)
        allocate(h(size(obj_y)))
        call har_negret_variance_path(obj_y, obj_x, coef, h)
        har_negret_nll = -gaussian_variance_loglik(obj_y(23:obj_ntrain), h(23:obj_ntrain)) / &
                         real(obj_ntrain - 22, dp)
        if (har_negret_nll /= har_negret_nll .or. har_negret_nll > 1.0e29_dp) har_negret_nll = 1.0e30_dp
        deallocate(h)
    end function har_negret_nll

    real(dp) function har_lev_nll(p)
        real(dp), intent(in) :: p(:)
        real(dp), allocatable :: h(:)
        real(dp) :: coef(7)

        coef = exp(p)
        allocate(h(size(obj_y)))
        call har_lev_variance_path(obj_y, obj_x, coef, h)
        har_lev_nll = -gaussian_variance_loglik(obj_y(23:obj_ntrain), h(23:obj_ntrain)) / &
                      real(obj_ntrain - 22, dp)
        if (har_lev_nll /= har_lev_nll .or. har_lev_nll > 1.0e29_dp) har_lev_nll = 1.0e30_dp
        deallocate(h)
    end function har_lev_nll

    real(dp) function semivar_har_nll(p)
        real(dp), intent(in) :: p(:)
        real(dp), allocatable :: h(:)
        real(dp) :: coef(7)

        coef = exp(p)
        allocate(h(size(obj_y)))
        call semivar_har_variance_path(obj_x, obj_x2, coef, h)
        semivar_har_nll = -gaussian_variance_loglik(obj_y(23:obj_ntrain), h(23:obj_ntrain)) / &
                          real(obj_ntrain - 22, dp)
        if (semivar_har_nll /= semivar_har_nll .or. semivar_har_nll > 1.0e29_dp) semivar_har_nll = 1.0e30_dp
        deallocate(h)
    end function semivar_har_nll

    real(dp) function midas_nll(p)
        real(dp), intent(in) :: p(:)
        real(dp), allocatable :: h(:)
        real(dp) :: a, bcoef, theta1, theta2

        a = exp(p(1))
        bcoef = exp(p(2))
        theta1 = exp(p(3))
        theta2 = exp(p(4))
        allocate(h(size(obj_y)))
        call midas_variance_path(obj_x, obj_k_lag, a, bcoef, theta1, theta2, h)
        midas_nll = -gaussian_variance_loglik(obj_y(obj_k_lag + 1:obj_ntrain), h(obj_k_lag + 1:obj_ntrain)) / &
                    real(obj_ntrain - obj_k_lag, dp)
        if (midas_nll /= midas_nll .or. midas_nll > 1.0e29_dp) midas_nll = 1.0e30_dp
        deallocate(h)
    end function midas_nll

    real(dp) function heavy_return_nll(p)
        real(dp), intent(in) :: p(:)
        real(dp), allocatable :: h(:)
        real(dp) :: omega, alpha, beta

        call heavy_transform(p, omega, alpha, beta)
        allocate(h(size(obj_y)))
        call heavy_variance_path(obj_x, omega, alpha, beta, obj_ntrain, h)
        heavy_return_nll = -gaussian_variance_loglik(obj_y(2:obj_ntrain), h(2:obj_ntrain)) / real(obj_ntrain - 1, dp)
        if (heavy_return_nll /= heavy_return_nll .or. heavy_return_nll > 1.0e29_dp) heavy_return_nll = 1.0e30_dp
        deallocate(h)
    end function heavy_return_nll

    real(dp) function heavy_rm_nll(p)
        real(dp), intent(in) :: p(:)
        real(dp) :: omega, alpha, beta

        call heavy_transform(p, omega, alpha, beta)
        heavy_rm_nll = -heavy_rm_loglik(omega, alpha, beta) / real(obj_ntrain - 1, dp)
        if (heavy_rm_nll /= heavy_rm_nll .or. heavy_rm_nll > 1.0e29_dp) heavy_rm_nll = 1.0e30_dp
    end function heavy_rm_nll

    real(dp) function ewma_affine_nll(p)
        real(dp), intent(in) :: p(:)
        real(dp), allocatable :: h(:)
        real(dp) :: lambda, a, b

        a = exp(p(1))
        b = exp(p(2))
        if (obj_fit_lambda) then
            lambda = lambda_transform(p(3))
        else
            lambda = obj_fixed_lambda
        end if
        allocate(h(size(obj_y)))
        call ewma_affine_variance_path(obj_x, lambda, a, b, obj_ntrain, h)
        ewma_affine_nll = -gaussian_variance_loglik(obj_y(1:obj_ntrain), h(1:obj_ntrain)) / real(obj_ntrain, dp)
        if (ewma_affine_nll /= ewma_affine_nll .or. ewma_affine_nll > 1.0e29_dp) ewma_affine_nll = 1.0e30_dp
        deallocate(h)
    end function ewma_affine_nll

    real(dp) function affine_nll(p)
        real(dp), intent(in) :: p(:)
        real(dp), allocatable :: h(:)
        real(dp) :: a, b

        a = exp(p(1))
        b = exp(p(2))
        allocate(h(size(obj_y)))
        call affine_variance_path(obj_x, a, b, h)
        affine_nll = -gaussian_variance_loglik(obj_y(1:obj_ntrain), h(1:obj_ntrain)) / real(obj_ntrain, dp)
        if (affine_nll /= affine_nll .or. affine_nll > 1.0e29_dp) affine_nll = 1.0e30_dp
        deallocate(h)
    end function affine_nll

    real(dp) function affine2_nll(p)
        real(dp), intent(in) :: p(:)
        real(dp), allocatable :: h(:)
        real(dp) :: a, b1, b2

        a = exp(p(1))
        b1 = exp(p(2))
        b2 = exp(p(3))
        allocate(h(size(obj_y)))
        call affine2_variance_path(obj_x, obj_x2, a, b1, b2, h)
        affine2_nll = -gaussian_variance_loglik(obj_y(1:obj_ntrain), h(1:obj_ntrain)) / real(obj_ntrain, dp)
        if (affine2_nll /= affine2_nll .or. affine2_nll > 1.0e29_dp) affine2_nll = 1.0e30_dp
        deallocate(h)
    end function affine2_nll

    real(dp) function heavy_rm_loglik(omega, alpha, beta)
        ! Gaussian quasi-likelihood for the HEAVY realized-measure equation.
        real(dp), intent(in) :: omega, alpha, beta
        real(dp), allocatable :: mu(:)
        integer :: t

        allocate(mu(size(obj_x)))
        call heavy_variance_path(obj_x, omega, alpha, beta, obj_ntrain, mu)
        heavy_rm_loglik = 0.0_dp
        do t = 2, obj_ntrain
            heavy_rm_loglik = heavy_rm_loglik - log_sqrt_2pi - 0.5_dp*log(max(mu(t), min_var)) - &
                              0.5_dp*max(obj_x(t), min_var) / max(mu(t), min_var)
        end do
        deallocate(mu)
    end function heavy_rm_loglik

    subroutine finite_diff_obj(nll_func, p, np, f, g)
        interface
            function nll_func(p) result(val)
                import :: dp
                real(dp), intent(in) :: p(:)
                real(dp) :: val
            end function nll_func
        end interface
        integer, intent(in) :: np
        real(dp), intent(in) :: p(np)
        real(dp), intent(out) :: f, g(np)
        real(dp), allocatable :: pp(:), pm(:)
        real(dp) :: step, fp, fm
        integer :: j

        f = nll_func(p)
        allocate(pp(np), pm(np))
        do j = 1, np
            step = 1.0e-5_dp * max(1.0_dp, abs(p(j)))
            pp = p
            pm = p
            pp(j) = pp(j) + step
            pm(j) = pm(j) - step
            fp = nll_func(pp)
            fm = nll_func(pm)
            g(j) = (fp - fm) / (2.0_dp*step)
        end do
        deallocate(pp, pm)
    end subroutine finite_diff_obj

    pure subroutine heavy_transform(p, omega, alpha, beta)
        ! Map unconstrained parameters to omega > 0, alpha >= 0, beta >= 0, alpha+beta < 0.999.
        real(dp), intent(in) :: p(:)
        real(dp), intent(out) :: omega, alpha, beta
        real(dp) :: persist, share

        omega = exp(p(1))
        persist = 0.999_dp / (1.0_dp + exp(-p(2)))
        share = 1.0_dp / (1.0_dp + exp(-p(3)))
        alpha = persist*share
        beta = persist*(1.0_dp - share)
    end subroutine heavy_transform

    pure subroutine heavy_inverse_transform(omega, alpha, beta, p)
        ! Map constrained HEAVY parameters to unconstrained optimizer parameters.
        real(dp), intent(in) :: omega, alpha, beta
        real(dp), intent(out) :: p(3)
        real(dp) :: persist, share

        persist = min(max(alpha + beta, 1.0e-8_dp), 0.998_dp)
        share = min(max(alpha / persist, 1.0e-8_dp), 1.0_dp - 1.0e-8_dp)
        p(1) = log(max(omega, min_var))
        p(2) = log((persist / 0.999_dp) / max(1.0_dp - persist / 0.999_dp, 1.0e-8_dp))
        p(3) = log(share / (1.0_dp - share))
    end subroutine heavy_inverse_transform

    pure real(dp) function lambda_transform(q)
        ! Unconstrained q to lambda in (0.5, 0.9999).
        real(dp), intent(in) :: q

        lambda_transform = 0.5_dp + 0.4999_dp / (1.0_dp + exp(-q))
    end function lambda_transform

    pure real(dp) function lambda_inverse(lambda)
        ! Lambda in (0.5,0.9999) to unconstrained q.
        real(dp), intent(in) :: lambda
        real(dp) :: x

        x = min(max((lambda - 0.5_dp) / 0.4999_dp, 1.0e-8_dp), 1.0_dp - 1.0e-8_dp)
        lambda_inverse = log(x / (1.0_dp - x))
    end function lambda_inverse

    pure subroutine har_row(x, t, row)
        ! Return [1, x_{t-1}, avg(x_{t-1..t-5}), avg(x_{t-1..t-22})].
        real(dp), intent(in) :: x(:)
        integer, intent(in) :: t
        real(dp), intent(out) :: row(4)

        row(1) = 1.0_dp
        row(2) = max(x(t - 1), min_var)
        row(3) = sum(max(x(t - 5:t - 1), min_var)) / 5.0_dp
        row(4) = sum(max(x(t - 22:t - 1), min_var)) / 22.0_dp
    end subroutine har_row

    pure subroutine harx_row(x, x_exog, t, row)
        ! Return HAR terms plus exogenous variance x_exog(t), known at forecast time.
        real(dp), intent(in) :: x(:), x_exog(:)
        integer, intent(in) :: t
        real(dp), intent(out) :: row(5)

        call har_row(x, t, row(1:4))
        row(5) = max(x_exog(t), min_var)
    end subroutine harx_row

    pure subroutine harx_lev_row(y, x, x_exog, t, row)
        ! Return HARX terms and daily/weekly/monthly negative-return-squared terms.
        real(dp), intent(in) :: y(:), x(:), x_exog(:)
        integer, intent(in) :: t
        real(dp), intent(out) :: row(8)

        call harx_row(x, x_exog, t, row(1:5))
        row(6) = min(y(t - 1), 0.0_dp)**2
        row(7) = sum(min(y(t - 5:t - 1), 0.0_dp)**2) / 5.0_dp
        row(8) = sum(min(y(t - 22:t - 1), 0.0_dp)**2) / 22.0_dp
    end subroutine harx_lev_row

    pure subroutine log_har_row(x, t, row)
        ! Return [1, log(x_{t-1}), log(avg_5(x)), log(avg_22(x))].
        real(dp), intent(in) :: x(:)
        integer, intent(in) :: t
        real(dp), intent(out) :: row(4)

        row(1) = 1.0_dp
        row(2) = log(max(x(t - 1), min_var))
        row(3) = log(max(sum(max(x(t - 5:t - 1), min_var)) / 5.0_dp, min_var))
        row(4) = log(max(sum(max(x(t - 22:t - 1), min_var)) / 22.0_dp, min_var))
    end subroutine log_har_row

    pure subroutine sqrt_har_row(x, t, row)
        ! Return [1, sqrt(x_{t-1}), sqrt(avg_5(x)), sqrt(avg_22(x))].
        real(dp), intent(in) :: x(:)
        integer, intent(in) :: t
        real(dp), intent(out) :: row(4)

        row(1) = 1.0_dp
        row(2) = sqrt(max(x(t - 1), min_var))
        row(3) = sqrt(max(sum(max(x(t - 5:t - 1), min_var)) / 5.0_dp, min_var))
        row(4) = sqrt(max(sum(max(x(t - 22:t - 1), min_var)) / 22.0_dp, min_var))
    end subroutine sqrt_har_row

    pure subroutine harj_row(x, jump, t, row)
        ! Return intercept, daily/weekly/monthly RV lags, and daily/weekly/monthly jump lags.
        real(dp), intent(in) :: x(:), jump(:)
        integer, intent(in) :: t
        real(dp), intent(out) :: row(7)

        row(1) = 1.0_dp
        row(2) = max(x(t - 1), min_var)
        row(3) = sum(max(x(t - 5:t - 1), min_var)) / 5.0_dp
        row(4) = sum(max(x(t - 22:t - 1), min_var)) / 22.0_dp
        row(5) = max(jump(t - 1), 0.0_dp)
        row(6) = sum(max(jump(t - 5:t - 1), 0.0_dp)) / 5.0_dp
        row(7) = sum(max(jump(t - 22:t - 1), 0.0_dp)) / 22.0_dp
    end subroutine harj_row

    pure subroutine har_lev_row(y, x, t, row)
        ! Return HAR terms and daily/weekly/monthly negative-return-squared leverage terms.
        real(dp), intent(in) :: y(:), x(:)
        integer, intent(in) :: t
        real(dp), intent(out) :: row(7)

        row(1) = 1.0_dp
        row(2) = max(x(t - 1), min_var)
        row(3) = sum(max(x(t - 5:t - 1), min_var)) / 5.0_dp
        row(4) = sum(max(x(t - 22:t - 1), min_var)) / 22.0_dp
        row(5) = min(y(t - 1), 0.0_dp)**2
        row(6) = sum(min(y(t - 5:t - 1), 0.0_dp)**2) / 5.0_dp
        row(7) = sum(min(y(t - 22:t - 1), 0.0_dp)**2) / 22.0_dp
    end subroutine har_lev_row

    pure subroutine semivar_har_row(x_pos, x_neg, t, row)
        ! Return intercept plus daily/weekly/monthly positive and negative semivariance lags.
        real(dp), intent(in) :: x_pos(:), x_neg(:)
        integer, intent(in) :: t
        real(dp), intent(out) :: row(7)

        row(1) = 1.0_dp
        row(2) = max(x_pos(t - 1), min_var)
        row(3) = sum(max(x_pos(t - 5:t - 1), min_var)) / 5.0_dp
        row(4) = sum(max(x_pos(t - 22:t - 1), min_var)) / 22.0_dp
        row(5) = max(x_neg(t - 1), min_var)
        row(6) = sum(max(x_neg(t - 5:t - 1), min_var)) / 5.0_dp
        row(7) = sum(max(x_neg(t - 22:t - 1), min_var)) / 22.0_dp
    end subroutine semivar_har_row

    pure subroutine beta_lag_weights(k_lag, theta1, theta2, weights)
        ! Normalized Beta-polynomial lag weights for lags 1..k_lag.
        integer, intent(in) :: k_lag
        real(dp), intent(in) :: theta1, theta2
        real(dp), intent(out) :: weights(k_lag)
        real(dp) :: x, logw, max_logw, sumw
        integer :: k

        max_logw = -huge(1.0_dp)
        do k = 1, k_lag
            x = real(k, dp) / real(k_lag + 1, dp)
            logw = (theta1 - 1.0_dp)*log(max(x, 1.0e-12_dp)) + &
                   (theta2 - 1.0_dp)*log(max(1.0_dp - x, 1.0e-12_dp))
            weights(k) = logw
            max_logw = max(max_logw, logw)
        end do
        sumw = 0.0_dp
        do k = 1, k_lag
            weights(k) = exp(weights(k) - max_logw)
            sumw = sumw + weights(k)
        end do
        if (sumw > 0.0_dp) then
            weights = weights / sumw
        else
            weights = 1.0_dp / real(k_lag, dp)
        end if
    end subroutine beta_lag_weights

end module realized_vol_forecast_mod

! GARCH model fitting and OOS evaluation for comparison with realized-variance models.
! Provides fit_garch_rv_model and eval_garch_oos_one_horizon.
!
! The GARCH conditional variance h_t (variance of daily log-returns) and the intraday
! realized variance RV_t can differ systematically in level AND in dynamic scaling.
! Overnight returns, bid-ask bounce, and illiquidity mean that for many assets the
! correlation between log(h_t) and log(RV_t) is well below 1, so the GARCH dynamics
! over-predict the variation in log(RV).
!
! Alignment: OLS regression of log(RV_t) on log(h_t) over the training sample:
!   log(RV_t) = ols_a + ols_b * log(h_t) + e_t
! Forecast: ols_a + ols_b * log(E[h_{t+k} | F_t])
! When ols_b << 1 the GARCH dynamics contribute little; when ols_b ~ 1 they matter fully.
!
! Multi-step forecast:
!   E[h_{t+k} | F_t] = omega/(1-p) + p^{k-1} * (h_{t+1|t} - omega/(1-p))
! where p is persistence and h_{t+1|t} = h(t+1) from the filtered variance path.
!
! Supported model names: "NAGARCH", "GJR", "GARCH", "REALGARCH".
! For REALGARCH the measurement equation E[log RV] = mu + phi*log h provides the direct
! log-RV mapping; no OLS alignment is applied.

module garch_rv_compare_mod
    use kind_mod,           only: dp
    use garch_types_mod,    only: garch_params_t
    use garch_fit_mod,      only: fit_nagarch,    nagarch_persist, &
                                   fit_gjr,         gjr_persist, &
                                   fit_symm_garch,  symm_garch_persist
    use garch_forecast_mod, only: nagarch_variance_path, gjr_variance_path, &
                                   symm_garch_variance_path
    use stats_mod,          only: simple_linreg
    use realized_garch_mod, only: fit_realized_garch, realized_garch_result_t
    implicit none
    private

    public :: garch_rv_model_t, fit_garch_rv_model, eval_garch_oos_one_horizon

    type :: garch_rv_model_t
        character(len=10)     :: name        = ""
        real(dp)              :: omega       = 0.0_dp
        real(dp)              :: persist     = 0.0_dp
        real(dp)              :: ols_a       = 0.0_dp  ! intercept: log(RV) = ols_a + ols_b*log(h)
        real(dp)              :: ols_b       = 0.0_dp  ! slope
        real(dp)              :: lh_inf     = 0.0_dp  ! fixed point of log-h recursion (REALGARCH only)
        real(dp), allocatable :: h(:)
    end type garch_rv_model_t

contains

    subroutine fit_garch_rv_model(ret, log_rv, ntrain, model_name, model)
        ! Fit a named GARCH model on training data and filter the full return series.
        ! ret(:)    daily close-to-close log-returns, length ndays
        ! log_rv(:) daily log realized variance (same length); used for OLS alignment
        real(dp),         intent(in)  :: ret(:)
        real(dp),         intent(in)  :: log_rv(:)
        integer,          intent(in)  :: ntrain
        character(len=*), intent(in)  :: model_name
        type(garch_rv_model_t), intent(out) :: model

        type(garch_params_t) :: params
        real(dp), allocatable :: log_h_train(:), x_rv(:)
        type(realized_garch_result_t) :: rg_result
        integer :: ndays
        real(dp) :: f_best
        integer :: niter
        logical :: conv
        integer, parameter :: max_iter = 1000
        real(dp), parameter :: gtol = 1.0e-6_dp

        ndays      = size(ret)
        model%name = model_name
        allocate(model%h(ndays))
        model%h    = 0.0_dp

        select case (trim(model_name))
        case ("NAGARCH")
            call fit_nagarch(ret(1:ntrain), max_iter, gtol, f_best, params, niter, conv)
            call nagarch_variance_path(ret, params, model%h)
            model%omega   = params%omega
            model%persist = nagarch_persist(params)
        case ("GJR")
            call fit_gjr(ret(1:ntrain), max_iter, gtol, f_best, params, niter, conv)
            call gjr_variance_path(ret, params, model%h)
            model%omega   = params%omega
            model%persist = gjr_persist(params)
        case ("GARCH")
            call fit_symm_garch(ret(1:ntrain), max_iter, gtol, f_best, params, niter, conv)
            call symm_garch_variance_path(ret, params, model%h)
            model%omega   = params%omega
            model%persist = symm_garch_persist(params)
        case ("REALGARCH")
            allocate(x_rv(ndays))
            x_rv = exp(log_rv)
            call fit_realized_garch(ret, x_rv, ntrain, max_iter, gtol, rg_result, model%h)
            deallocate(x_rv)
            model%omega  = rg_result%params%omega
            model%persist = rg_result%persist
            model%ols_a  = rg_result%params%mu   ! measurement eq: E[log RV] = mu + phi*log h
            model%ols_b  = rg_result%params%phi
            model%lh_inf = (rg_result%params%omega + rg_result%params%beta * rg_result%params%mu) / &
                            max(1.0_dp - rg_result%persist, 1.0e-8_dp)
            return  ! measurement equation provides exact log-RV mapping; no OLS alignment needed
        case default
            write(*,'(2A)') "fit_garch_rv_model: unknown model ", trim(model_name)
            return
        end select

        ! OLS alignment on training sample: log(RV_t) = ols_a + ols_b * log(h_t)
        ! ols_b << 1 means GARCH dynamics are weakly correlated with log(RV) dynamics
        allocate(log_h_train(ntrain))
        log_h_train = log(max(model%h(1:ntrain), 1.0e-20_dp))
        call simple_linreg(log_h_train, log_rv(1:ntrain), model%ols_b, model%ols_a)
        deallocate(log_h_train)
    end subroutine fit_garch_rv_model


    subroutine eval_garch_oos_one_horizon(log_rv, ndays, ntrain, h_step, model, rmse_out, qlike_out)
        ! OOS RMSE and exp-QLIKE for a fitted GARCH model at forecast horizon h_step.
        ! Log-variance forecast: ols_a + ols_b * log(E[h_{t+h} | F_t])
        real(dp), intent(in) :: log_rv(:)
        integer,  intent(in) :: ndays, ntrain, h_step
        type(garch_rv_model_t), intent(in)  :: model
        real(dp), intent(out) :: rmse_out, qlike_out

        integer  :: t, neval
        real(dp) :: actual, fc_logvar, e, var_inf, h_fc, lh_1, lh_k

        rmse_out  = 0.0_dp
        qlike_out = 0.0_dp
        neval     = 0

        if (.not. allocated(model%h) .or. model%persist <= 0.0_dp) then
            rmse_out  = huge(1.0_dp)
            qlike_out = huge(1.0_dp)
            return
        end if

        var_inf = model%omega / max(1.0_dp - model%persist, 1.0e-8_dp)

        do t = ntrain + 1, ndays - h_step
            actual = log_rv(t + h_step)
            if (model%name == "REALGARCH") then
                ! log-h is AR(1): E[log h_{t+k}] = lh_inf + persist^{k-1}*(log h_{t+1} - lh_inf)
                ! Measurement eq: E[log RV_{t+k}] = mu + phi * E[log h_{t+k}]
                lh_1      = log(max(model%h(t + 1), 1.0e-20_dp))
                lh_k      = model%lh_inf + model%persist**(h_step - 1) * (lh_1 - model%lh_inf)
                fc_logvar = model%ols_a + model%ols_b * lh_k
            else
                ! h level AR(1): E[h_{t+k}] = var_inf + persist^{k-1}*(h_{t+1} - var_inf)
                h_fc      = var_inf + model%persist**(h_step - 1) * (model%h(t + 1) - var_inf)
                fc_logvar = model%ols_a + model%ols_b * log(max(h_fc, 1.0e-20_dp))
            end if
            e         = actual - fc_logvar
            rmse_out  = rmse_out  + e**2
            qlike_out = qlike_out + (exp(e) - e - 1.0_dp)
            neval     = neval + 1
        end do

        if (neval > 0) then
            rmse_out  = sqrt(rmse_out  / real(neval, dp))
            qlike_out = qlike_out / real(neval, dp)
        else
            rmse_out  = huge(1.0_dp)
            qlike_out = huge(1.0_dp)
        end if
    end subroutine eval_garch_oos_one_horizon

end module garch_rv_compare_mod

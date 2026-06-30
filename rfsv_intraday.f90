! RFSV intraday diagnostics and OOS forecast evaluation routines.

module rfsv_intraday_mod
    use kind_mod,                  only: dp
    use stats_mod,                 only: fit_ar_direct
    use realized_vol_forecast_mod, only: fit_log_har_direct, fit_log_har_asym_direct
    use rough_sv_mod,              only: rfsv_forecast_weights
    implicit none
    private

    public :: print_log_msd_table, eval_oos_one_horizon

contains

    subroutine print_log_msd_table(x, lags, h_hat, sigma_hat)
        ! Print observed vs fitted log-MSD table for RFSV model diagnostics.
        real(dp), intent(in) :: x(:)        ! daily log realized variance series
        integer,  intent(in) :: lags(:)     ! lags (in days) at which to evaluate log-MSD
        real(dp), intent(in) :: h_hat       ! estimated Hurst roughness index
        real(dp), intent(in) :: sigma_hat   ! estimated vol of log-RV

        integer  :: i, k, n
        real(dp) :: msd, log_msd_data, log_msd_fit

        n = size(x)
        write(*,'(A)') "Log-MSD diagnostic (slope = 2H):"
        write(*,'(A5,2(1X,A12),1X,A8)') "lag", "log-MSD_data", "log-MSD_fit", "ratio"
        write(*,'(A)') repeat("-", 44)
        do i = 1, size(lags)
            k = lags(i)
            if (k >= n) cycle
            msd = sum((x(k+1:n) - x(1:n-k))**2) / real(n-k, dp)
            if (msd <= 0.0_dp) cycle
            log_msd_data = log(msd)
            log_msd_fit  = log(sigma_hat**2) + 2.0_dp * h_hat * log(real(k, dp))
            write(*,'(I5,2(1X,F12.4),1X,F8.4)') &
                k, log_msd_data, log_msd_fit, exp(log_msd_data - log_msd_fit)
        end do
        write(*,'(A)') repeat("-", 44)
        write(*,*)
    end subroutine print_log_msd_table

    subroutine eval_oos_one_horizon(log_rv, log_rsv_neg, ndays, ntrain, p, h, use_rfsv, hurst, &
                                     rmse_ar_out, qlike_ar_out, &
                                     rmse_har_out, qlike_har_out, &
                                     rmse_hara_out, qlike_hara_out, &
                                     rmse_rfsv_out, qlike_rfsv_out, &
                                     print_row)
        ! Evaluate AR(p), log-HAR, HAR-A, and optionally RFSV OOS forecasts for horizon h.
        ! Prints one result row per model column set unless print_row=.false.
        real(dp), intent(in)            :: log_rv(:)       ! full daily log realized variance series
        real(dp), intent(in)            :: log_rsv_neg(:)  ! full daily log negative semi-variance
        integer,  intent(in)            :: ndays           ! total number of days
        integer,  intent(in)            :: ntrain          ! number of training days
        integer,  intent(in)            :: p               ! AR order (also used as RFSV lookback)
        integer,  intent(in)            :: h               ! forecast horizon (days)
        logical,  intent(in)            :: use_rfsv        ! if true, compute RFSV model forecast
        real(dp), intent(in)            :: hurst           ! estimated Hurst index (used when use_rfsv)
        real(dp), intent(out), optional :: rmse_ar_out     ! RMSE for AR(p) model
        real(dp), intent(out), optional :: qlike_ar_out    ! QLIKE for AR(p) model
        real(dp), intent(out), optional :: rmse_har_out    ! RMSE for log-HAR model
        real(dp), intent(out), optional :: qlike_har_out   ! QLIKE for log-HAR model
        real(dp), intent(out), optional :: rmse_hara_out   ! RMSE for HAR-A model
        real(dp), intent(out), optional :: qlike_hara_out  ! QLIKE for HAR-A model
        real(dp), intent(out), optional :: rmse_rfsv_out   ! RMSE for RFSV model
        real(dp), intent(out), optional :: qlike_rfsv_out  ! QLIKE for RFSV model
        logical,  intent(in), optional  :: print_row       ! if .false., suppress the output row

        real(dp) :: nu, beta_ar(p), beta_har(3), beta_hara(4), fc_ar, fc_har, fc_hara
        real(dp) :: rmse_ar, rmse_har, rmse_hara, qlike_ar, qlike_har, qlike_hara
        real(dp) :: e_ar, e_har, e_hara, actual
        integer  :: t, i, nfit_ar, nfit_har, nfit_hara, neval
        logical  :: do_print

        real(dp), allocatable :: beta_rfsv(:), delta_x(:)
        real(dp) :: fc_rfsv, rmse_rfsv, qlike_rfsv, e_rfsv

        do_print = .true.
        if (present(print_row)) do_print = print_row

        nu = sum(log_rv(1:ntrain)) / real(ntrain, dp)
        call fit_ar_direct(log_rv(1:ntrain) - nu, p, h, beta_ar, nfit_ar)
        call fit_log_har_direct(log_rv(1:ntrain), h, beta_har, nfit_har)
        call fit_log_har_asym_direct(log_rv(1:ntrain), log_rsv_neg(1:ntrain), h, beta_hara, nfit_hara)

        rmse_ar   = 0.0_dp; qlike_ar   = 0.0_dp
        rmse_har  = 0.0_dp; qlike_har  = 0.0_dp
        rmse_hara = 0.0_dp; qlike_hara = 0.0_dp
        rmse_rfsv = 0.0_dp; qlike_rfsv = 0.0_dp
        neval = 0

        allocate(beta_rfsv(p), delta_x(p))
        beta_rfsv = 0.0_dp; delta_x = 0.0_dp
        if (use_rfsv) call rfsv_forecast_weights(h, p, hurst, beta_rfsv)

        do t = ntrain + 1, ndays - h
            actual   = log_rv(t + h)
            fc_ar    = nu + sum(beta_ar * (log_rv(t:t-p+1:-1) - nu))
            fc_har   = beta_har(1) &
                     + beta_har(2) * log_rv(t) &
                     + beta_har(3) * sum(log_rv(max(1,t-21):t)) / real(min(22, t), dp)
            fc_hara  = beta_hara(1) &
                     + beta_hara(2) * log_rv(t) &
                     + beta_hara(3) * log_rsv_neg(t) &
                     + beta_hara(4) * sum(log_rv(max(1,t-21):t)) / real(min(22, t), dp)
            neval    = neval + 1
            e_ar     = actual - fc_ar
            e_har    = actual - fc_har
            e_hara   = actual - fc_hara
            rmse_ar   = rmse_ar   + e_ar**2
            rmse_har  = rmse_har  + e_har**2
            rmse_hara = rmse_hara + e_hara**2
            qlike_ar   = qlike_ar   + (exp(e_ar)   - e_ar   - 1.0_dp)
            qlike_har  = qlike_har  + (exp(e_har)  - e_har  - 1.0_dp)
            qlike_hara = qlike_hara + (exp(e_hara) - e_hara - 1.0_dp)

            if (use_rfsv) then
                do i = 1, p
                    delta_x(i) = log_rv(t-i+1) - log_rv(t-i)
                end do
                fc_rfsv    = log_rv(t) + dot_product(beta_rfsv, delta_x)
                e_rfsv     = actual - fc_rfsv
                rmse_rfsv  = rmse_rfsv  + e_rfsv**2
                qlike_rfsv = qlike_rfsv + (exp(e_rfsv) - e_rfsv - 1.0_dp)
            end if
        end do

        if (neval > 0) then
            rmse_ar    = sqrt(rmse_ar    / real(neval, dp))
            rmse_har   = sqrt(rmse_har   / real(neval, dp))
            rmse_hara  = sqrt(rmse_hara  / real(neval, dp))
            qlike_ar   = qlike_ar   / real(neval, dp)
            qlike_har  = qlike_har  / real(neval, dp)
            qlike_hara = qlike_hara / real(neval, dp)
            if (use_rfsv) then
                rmse_rfsv  = sqrt(rmse_rfsv  / real(neval, dp))
                qlike_rfsv = qlike_rfsv / real(neval, dp)
            end if
        end if

        if (do_print) then
            if (use_rfsv) then
                write(*,'(I6,4(2(1X,F10.4)))') h, &
                    rmse_ar, qlike_ar, rmse_har, qlike_har, rmse_hara, qlike_hara, rmse_rfsv, qlike_rfsv
            else
                write(*,'(I6,3(2(1X,F10.4)))') h, &
                    rmse_ar, qlike_ar, rmse_har, qlike_har, rmse_hara, qlike_hara
            end if
        end if

        if (present(rmse_ar_out))    rmse_ar_out    = rmse_ar
        if (present(qlike_ar_out))   qlike_ar_out   = qlike_ar
        if (present(rmse_har_out))   rmse_har_out   = rmse_har
        if (present(qlike_har_out))  qlike_har_out  = qlike_har
        if (present(rmse_hara_out))  rmse_hara_out  = rmse_hara
        if (present(qlike_hara_out)) qlike_hara_out = qlike_hara
        if (present(rmse_rfsv_out))  rmse_rfsv_out  = rmse_rfsv
        if (present(qlike_rfsv_out)) qlike_rfsv_out = qlike_rfsv
        deallocate(beta_rfsv, delta_x)
    end subroutine eval_oos_one_horizon

end module rfsv_intraday_mod

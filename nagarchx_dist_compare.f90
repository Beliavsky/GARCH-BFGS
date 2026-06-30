! Driver module for comparing NAGARCH and NAGARCH-X models across noise distributions.
!
! Models compared per asset:
!   NAGARCH (Normal)     -- k=4 baseline, always fitted
!   NAGARCH-X <dist>     -- k depends on distribution, fitted for each entry in dist_names
!
! dist_names is passed in by the caller (typically xfit_nagarchx_dist.f90).
! Supported distributions: "NORMAL" (k=5), "T" (k=6), "FS_SKEWT" (k=7).
! Warm-start chain (when prerequisites are present in dist_names):
!   NAGARCH -> NAGARCH-X(N) -> NAGARCH-X(T) -> NAGARCH-X(S)
!
! Two analyses per asset:
!   Train/test split  -- logL, QLIKE, QLIKE_RV, Corr(h,RV) on held-out test set
!   Full sample       -- logL, AIC, BIC on all observations

module nagarchx_dist_compare_mod
    use kind_mod,        only: dp
    use date_mod,        only: date_label
    use market_data_mod, only: ohlcv_series_t, read_intraday_prices_file, filter_intraday_session
    use intraday_realized_measures_mod, only: daily_realized_panel_t, build_daily_realized_panel
    use garch_types_mod,    only: garch_params_t
    use garch_fit_mod,      only: fit_nagarch, nagarch_persist
    use garch_forecast_mod, only: nagarch_variance_path
    use stats_mod,          only: mean, demean_first, gaussian_loglik, correlation
    use program_utils_mod,  only: elapsed_since
    use path_utils_mod,     only: asset_label
    use distributions_mod,  only: dist_t, logpdf_std, pdf_fs_skewt
    use nagarchx_mod,       only: nagarchx_params_t, fit_nagarchx, nagarchx_variance_path, nagarchx_persist
    use nagarchx_t_mod,     only: nagarchx_t_params_t, fit_nagarchx_t, nagarchx_t_persist
    use nagarchx_skewt_mod, only: nagarchx_skewt_params_t, fit_nagarchx_skewt, nagarchx_skewt_persist
    implicit none
    private

    real(dp), parameter :: min_var       = 1.0e-12_dp
    real(dp), parameter :: min_pdf       = 1.0e-300_dp
    integer,  parameter :: default_ntest = 250
    integer,  parameter :: max_iter      = 500
    real(dp), parameter :: gtol          = 1.0e-5_dp
    real(dp), parameter :: trading_days  = 252.0_dp
    integer,  parameter :: k_nag         = 4

    ! Per-asset result; dist-indexed arrays are allocated to size(dist_names) by compare_one_file.
    type, public :: nagarchx_dist_result_t
        character(len=8) :: asset  = ""
        integer          :: ntrain = 0
        integer          :: ntest  = 0
        integer          :: nobs   = 0
        ! NAGARCH (Normal) baseline -- train/test split
        real(dp) :: nag_logL_tr  = 0.0_dp
        real(dp) :: nag_logL_te  = 0.0_dp
        real(dp) :: nag_aic      = 0.0_dp
        real(dp) :: nag_bic      = 0.0_dp
        real(dp) :: nag_qlike    = 0.0_dp
        real(dp) :: nag_qlike_rv = 0.0_dp
        real(dp) :: nag_corr_rv  = 0.0_dp
        real(dp) :: nag_vol      = 0.0_dp
        integer  :: nag_niter    = 0
        logical  :: nag_conv     = .false.
        ! NAGARCH (Normal) baseline -- full sample
        real(dp) :: nag_logL_full = 0.0_dp
        real(dp) :: nag_aic_full  = 0.0_dp
        real(dp) :: nag_bic_full  = 0.0_dp
        ! NAGARCH-X per distribution -- train/test split (size = size(dist_names))
        real(dp), allocatable :: logL_tr(:), logL_te(:), aic(:), bic(:)
        real(dp), allocatable :: qlike(:), qlike_rv(:), corr_rv(:), vol(:)
        real(dp), allocatable :: nu(:), xi(:)
        integer,  allocatable :: niter(:)
        logical,  allocatable :: conv(:)
        ! NAGARCH-X per distribution -- full sample
        real(dp), allocatable :: logL_full(:), aic_full(:), bic_full(:)
    end type nagarchx_dist_result_t

    public :: run_nagarchx_dist_compare

contains

    subroutine run_nagarchx_dist_compare(filenames, dist_names)
        character(len=*), intent(in) :: filenames(:)
        character(len=*), intent(in) :: dist_names(:)
        type(nagarchx_dist_result_t), allocatable :: results(:)
        integer :: i, id

        write(*, '(A)', advance='no') "Distribution abbreviations:"
        do id = 1, size(dist_names)
            write(*, '(2X,A,A,A)', advance='no') &
                  trim(dist_abbrev(dist_names(id))), "=", trim(dist_fullname(dist_names(id)))
        end do
        write(*, '()')

        allocate(results(size(filenames)))
        do i = 1, size(filenames)
            if (i > 1) print '(A)', ""
            call compare_one_file(trim(filenames(i)), dist_names, results(i))
        end do
        if (size(filenames) > 1) call print_asset_summary(results, dist_names)
        deallocate(results)
    end subroutine run_nagarchx_dist_compare

    subroutine compare_one_file(filename, dist_names, result)
        character(len=*),             intent(in)  :: filename
        character(len=*),             intent(in)  :: dist_names(:)
        type(nagarchx_dist_result_t), intent(out) :: result
        type(ohlcv_series_t) :: bars, regular_bars
        type(daily_realized_panel_t) :: daily
        real(dp), allocatable :: ret_cc(:), rv(:), h_nag(:), h_nagx(:,:)
        real(dp), allocatable :: h_nag_fs(:), h_nagx_fs(:,:)
        type(garch_params_t)          :: p_nag, p_nag_fs
        type(nagarchx_params_t)       :: p_nagx_n, p_warm_n
        type(nagarchx_t_params_t)     :: p_nagx_t, p_warm_t, p_warm_t_fs
        type(nagarchx_skewt_params_t) :: p_nagx_s, p_warm_s, p_warm_s_fs
        type(nagarchx_params_t)       :: p_nagx_n_fs
        type(nagarchx_t_params_t)     :: p_nagx_t_fs
        type(nagarchx_skewt_params_t) :: p_nagx_s_fs
        real(dp) :: fopt, persist, mu_y2, mu_rv, logn, logn_n
        real(dp) :: t0, t1, read_sec, fit_sec, elapsed_sec
        integer  :: nobs, ntrain, test_first, te_last, ndist, id, k
        integer  :: niter_tmp
        logical  :: conv_tmp
        logical  :: fitted_normal, fitted_t, fitted_n_fs, fitted_t_fs

        ndist = size(dist_names)

        call cpu_time(t0)
        call cpu_time(t1)
        call read_intraday_prices_file(filename, bars)
        call filter_intraday_session(bars, regular_bars)
        read_sec = elapsed_since(t1)

        call build_daily_realized_panel(regular_bars, daily)
        nobs = size(daily%date) - 1
        if (nobs < default_ntest + 30) error stop "compare_one_file: not enough daily observations"
        ntrain     = nobs - default_ntest
        test_first = ntrain + 1
        te_last    = nobs

        allocate(ret_cc(nobs), rv(nobs), h_nag(nobs), h_nagx(nobs, ndist))
        allocate(h_nag_fs(nobs), h_nagx_fs(nobs, ndist))
        ret_cc = log(daily%close(2:nobs+1) / daily%close(1:nobs))
        rv     = daily%rv(1:nobs)
        call demean_first(ret_cc, ntrain)

        print '(A,A)',          "Input file:   ", trim(filename)
        print '(A,I0,A,A,A,A)', "Daily obs:    ", nobs, &
              " from ", date_label(daily%date(1)), " to ", date_label(daily%date(nobs+1))
        print '(A,I0,A,A,A,A)', "Training:     ", ntrain, &
              " through ", date_label(daily%date(ntrain+1)), &
              ";  test starts ", date_label(daily%date(ntrain+2))
        print '(A,I0)',          "Test obs:     ", default_ntest
        print '(A,ES10.3,A,ES10.3)', "Mean(y^2):    ", mean(ret_cc(1:ntrain)**2), &
              "  Mean(RV): ", mean(rv(1:ntrain))

        mu_y2 = max(mean(ret_cc(1:ntrain)**2), min_var)
        mu_rv  = max(mean(rv(1:ntrain)),         min_var)
        logn   = log(real(ntrain, dp))
        logn_n = log(real(nobs,   dp))

        ! Allocate result arrays
        allocate(result%logL_tr(ndist), result%logL_te(ndist), result%aic(ndist), result%bic(ndist))
        allocate(result%qlike(ndist), result%qlike_rv(ndist), result%corr_rv(ndist), result%vol(ndist))
        allocate(result%nu(ndist), result%xi(ndist), result%niter(ndist), result%conv(ndist))
        allocate(result%logL_full(ndist), result%aic_full(ndist), result%bic_full(ndist))
        result%nu = 0.0_dp;  result%xi = 0.0_dp

        call cpu_time(t1)

        !--------------------------------------------------------------------
        ! NAGARCH(Normal) baseline -- train
        !--------------------------------------------------------------------
        call fit_nagarch(ret_cc(1:ntrain), max_iter, gtol, fopt, p_nag, &
                         result%nag_niter, result%nag_conv)
        call nagarch_variance_path(ret_cc, p_nag, h_nag)
        persist = nagarch_persist(p_nag)

        print '(/,A)', "NAGARCH parameters:"
        print '(3X,A,ES10.3,A,F8.4,A,F8.4,A,F7.3,A,F7.4,A,I0,A,L1)', &
              "omega=", p_nag%omega, "  alpha=", p_nag%alpha, &
              "  beta=",  p_nag%beta,  "  theta=", p_nag%theta, &
              "  persist=", persist, "  iter=", result%nag_niter, "  conv=", result%nag_conv

        !--------------------------------------------------------------------
        ! NAGARCH-X -- train, one distribution at a time
        !--------------------------------------------------------------------
        fitted_normal = .false.
        fitted_t      = .false.

        do id = 1, ndist
            select case (trim(dist_names(id)))

            case ("NORMAL")
                p_warm_n = nagarchx_params_t(omega = p_nag%omega, alpha = p_nag%alpha, &
                                             beta  = p_nag%beta,  theta = p_nag%theta, &
                                             delta = 0.05_dp * mu_y2 / mu_rv)
                call fit_nagarchx(ret_cc, rv, ntrain, max_iter, gtol, fopt, p_nagx_n, &
                                  result%niter(id), result%conv(id), warm=p_warm_n)
                call nagarchx_variance_path(ret_cc, rv, p_nagx_n, h_nagx(:,id))
                persist = nagarchx_persist(p_nagx_n)
                fitted_normal = .true.
                print '(/,A)', "NAGARCH-X N parameters:"
                print '(3X,A,ES10.3,A,F8.4,A,F8.4,A,F7.3,A,ES10.3,A,F7.4,A,I0,A,L1)', &
                      "omega=", p_nagx_n%omega, "  alpha=", p_nagx_n%alpha, &
                      "  beta=",  p_nagx_n%beta,  "  theta=", p_nagx_n%theta, &
                      "  delta=", p_nagx_n%delta, "  persist=", persist, &
                      "  iter=", result%niter(id), "  conv=", result%conv(id)

            case ("T")
                if (fitted_normal) then
                    p_warm_t = nagarchx_t_params_t(garch=p_nagx_n, nu=8.0_dp)
                    call fit_nagarchx_t(ret_cc, rv, ntrain, max_iter, gtol, fopt, p_nagx_t, &
                                        result%niter(id), result%conv(id), warm=p_warm_t)
                else
                    call fit_nagarchx_t(ret_cc, rv, ntrain, max_iter, gtol, fopt, p_nagx_t, &
                                        result%niter(id), result%conv(id))
                end if
                call nagarchx_variance_path(ret_cc, rv, p_nagx_t%garch, h_nagx(:,id))
                persist = nagarchx_t_persist(p_nagx_t)
                result%nu(id) = p_nagx_t%nu
                fitted_t = .true.
                print '(/,A)', "NAGARCH-X T parameters:"
                print '(3X,A,ES10.3,A,F8.4,A,F8.4,A,F7.3,A,ES10.3,A,F6.2,A,F7.4,A,I0,A,L1)', &
                      "omega=", p_nagx_t%garch%omega, "  alpha=", p_nagx_t%garch%alpha, &
                      "  beta=",  p_nagx_t%garch%beta,  "  theta=", p_nagx_t%garch%theta, &
                      "  delta=", p_nagx_t%garch%delta, "  nu=", p_nagx_t%nu, &
                      "  persist=", persist, "  iter=", result%niter(id), "  conv=", result%conv(id)

            case ("FS_SKEWT")
                if (fitted_t) then
                    p_warm_s = nagarchx_skewt_params_t(garch=p_nagx_t%garch, nu=p_nagx_t%nu, xi=1.0_dp)
                    call fit_nagarchx_skewt(ret_cc, rv, ntrain, max_iter, gtol, fopt, p_nagx_s, &
                                            result%niter(id), result%conv(id), warm=p_warm_s)
                else if (fitted_normal) then
                    p_warm_s = nagarchx_skewt_params_t(garch=p_nagx_n, nu=8.0_dp, xi=1.0_dp)
                    call fit_nagarchx_skewt(ret_cc, rv, ntrain, max_iter, gtol, fopt, p_nagx_s, &
                                            result%niter(id), result%conv(id), warm=p_warm_s)
                else
                    call fit_nagarchx_skewt(ret_cc, rv, ntrain, max_iter, gtol, fopt, p_nagx_s, &
                                            result%niter(id), result%conv(id))
                end if
                call nagarchx_variance_path(ret_cc, rv, p_nagx_s%garch, h_nagx(:,id))
                persist = nagarchx_skewt_persist(p_nagx_s)
                result%nu(id) = p_nagx_s%nu
                result%xi(id) = p_nagx_s%xi
                print '(/,A)', "NAGARCH-X S parameters:"
                print '(3X,A,ES10.3,A,F8.4,A,F8.4,A,F7.3,A,ES10.3,A,F6.2,A,F8.4,A,F7.4,A,I0,A,L1)', &
                      "omega=", p_nagx_s%garch%omega, "  alpha=", p_nagx_s%garch%alpha, &
                      "  beta=",  p_nagx_s%garch%beta,  "  theta=", p_nagx_s%garch%theta, &
                      "  delta=", p_nagx_s%garch%delta, "  nu=", p_nagx_s%nu, &
                      "  xi=", p_nagx_s%xi, "  persist=", persist, &
                      "  iter=", result%niter(id), "  conv=", result%conv(id)

            case default
                write(*, '(A,A)') "Warning: unknown distribution ignored: ", trim(dist_names(id))
                h_nagx(:,id) = h_nag

            end select
        end do

        fit_sec     = elapsed_since(t1)
        elapsed_sec = elapsed_since(t0)

        !--------------------------------------------------------------------
        ! Train/test metrics
        !--------------------------------------------------------------------
        result%asset  = asset_label(filename)
        result%ntrain = ntrain
        result%ntest  = default_ntest
        result%nobs   = nobs

        result%nag_logL_tr  = gaussian_loglik(ret_cc(1:ntrain), h_nag(1:ntrain))
        result%nag_logL_te  = gaussian_loglik(ret_cc(test_first:te_last), h_nag(test_first:te_last))
        result%nag_aic      = 2.0_dp*k_nag - 2.0_dp*result%nag_logL_tr
        result%nag_bic      = logn*k_nag   - 2.0_dp*result%nag_logL_tr
        result%nag_qlike    = mean(log(max(h_nag(test_first:te_last), min_var)) + &
                                   ret_cc(test_first:te_last)**2 / max(h_nag(test_first:te_last), min_var))
        result%nag_qlike_rv = mean(log(max(h_nag(test_first:te_last), min_var)) + &
                                   rv(test_first:te_last) / max(h_nag(test_first:te_last), min_var))
        result%nag_corr_rv  = correlation(h_nag(test_first:te_last), rv(test_first:te_last))
        result%nag_vol      = 100.0_dp * sqrt(trading_days * mean(h_nag(test_first:te_last)))

        do id = 1, ndist
            k = dist_nparams(dist_names(id))
            select case (trim(dist_names(id)))
            case ("NORMAL")
                result%logL_tr(id) = gaussian_loglik(ret_cc(1:ntrain), h_nagx(1:ntrain,id))
                result%logL_te(id) = gaussian_loglik(ret_cc(test_first:te_last), h_nagx(test_first:te_last,id))
            case ("T")
                result%logL_tr(id) = t_loglik(ret_cc(1:ntrain), h_nagx(1:ntrain,id), result%nu(id))
                result%logL_te(id) = t_loglik(ret_cc(test_first:te_last), h_nagx(test_first:te_last,id), result%nu(id))
            case ("FS_SKEWT")
                result%logL_tr(id) = skewt_loglik(ret_cc(1:ntrain), h_nagx(1:ntrain,id), result%nu(id), result%xi(id))
                result%logL_te(id) = skewt_loglik(ret_cc(test_first:te_last), h_nagx(test_first:te_last,id), &
                                                   result%nu(id), result%xi(id))
            case default
                result%logL_tr(id) = result%nag_logL_tr
                result%logL_te(id) = result%nag_logL_te
            end select
            result%aic(id)      = 2.0_dp*k - 2.0_dp*result%logL_tr(id)
            result%bic(id)      = logn*k    - 2.0_dp*result%logL_tr(id)
            result%qlike(id)    = mean(log(max(h_nagx(test_first:te_last,id), min_var)) + &
                                       ret_cc(test_first:te_last)**2 / max(h_nagx(test_first:te_last,id), min_var))
            result%qlike_rv(id) = mean(log(max(h_nagx(test_first:te_last,id), min_var)) + &
                                       rv(test_first:te_last) / max(h_nagx(test_first:te_last,id), min_var))
            result%corr_rv(id)  = correlation(h_nagx(test_first:te_last,id), rv(test_first:te_last))
            result%vol(id)      = 100.0_dp * sqrt(trading_days * mean(h_nagx(test_first:te_last,id)))
        end do

        !--------------------------------------------------------------------
        ! Full-sample fits, warm-started from training-set results
        !--------------------------------------------------------------------
        call fit_nagarch(ret_cc(1:nobs), max_iter, gtol, fopt, p_nag_fs, niter_tmp, conv_tmp)
        call nagarch_variance_path(ret_cc, p_nag_fs, h_nag_fs)
        result%nag_logL_full = gaussian_loglik(ret_cc(1:nobs), h_nag_fs(1:nobs))
        result%nag_aic_full  = 2.0_dp*k_nag - 2.0_dp*result%nag_logL_full
        result%nag_bic_full  = logn_n*k_nag  - 2.0_dp*result%nag_logL_full

        fitted_n_fs = .false.
        fitted_t_fs = .false.

        do id = 1, ndist
            k = dist_nparams(dist_names(id))
            select case (trim(dist_names(id)))

            case ("NORMAL")
                if (fitted_normal) then
                    call fit_nagarchx(ret_cc, rv, nobs, max_iter, gtol, fopt, p_nagx_n_fs, &
                                      niter_tmp, conv_tmp, warm=p_nagx_n)
                else
                    call fit_nagarchx(ret_cc, rv, nobs, max_iter, gtol, fopt, p_nagx_n_fs, &
                                      niter_tmp, conv_tmp)
                end if
                call nagarchx_variance_path(ret_cc, rv, p_nagx_n_fs, h_nagx_fs(:,id))
                result%logL_full(id) = gaussian_loglik(ret_cc(1:nobs), h_nagx_fs(1:nobs,id))
                fitted_n_fs = .true.

            case ("T")
                if (fitted_n_fs) then
                    p_warm_t_fs = nagarchx_t_params_t(garch=p_nagx_n_fs, &
                                                       nu=merge(p_nagx_t%nu, 8.0_dp, fitted_t))
                    call fit_nagarchx_t(ret_cc, rv, nobs, max_iter, gtol, fopt, p_nagx_t_fs, &
                                        niter_tmp, conv_tmp, warm=p_warm_t_fs)
                else if (fitted_t) then
                    call fit_nagarchx_t(ret_cc, rv, nobs, max_iter, gtol, fopt, p_nagx_t_fs, &
                                        niter_tmp, conv_tmp, warm=p_nagx_t)
                else
                    call fit_nagarchx_t(ret_cc, rv, nobs, max_iter, gtol, fopt, p_nagx_t_fs, &
                                        niter_tmp, conv_tmp)
                end if
                call nagarchx_variance_path(ret_cc, rv, p_nagx_t_fs%garch, h_nagx_fs(:,id))
                result%logL_full(id) = t_loglik(ret_cc(1:nobs), h_nagx_fs(1:nobs,id), p_nagx_t_fs%nu)
                fitted_t_fs = .true.

            case ("FS_SKEWT")
                if (fitted_t_fs) then
                    p_warm_s_fs = nagarchx_skewt_params_t(garch=p_nagx_t_fs%garch, &
                                                           nu=p_nagx_t_fs%nu, xi=1.0_dp)
                    call fit_nagarchx_skewt(ret_cc, rv, nobs, max_iter, gtol, fopt, p_nagx_s_fs, &
                                            niter_tmp, conv_tmp, warm=p_warm_s_fs)
                else if (fitted_n_fs) then
                    p_warm_s_fs = nagarchx_skewt_params_t(garch=p_nagx_n_fs, nu=8.0_dp, xi=1.0_dp)
                    call fit_nagarchx_skewt(ret_cc, rv, nobs, max_iter, gtol, fopt, p_nagx_s_fs, &
                                            niter_tmp, conv_tmp, warm=p_warm_s_fs)
                else
                    call fit_nagarchx_skewt(ret_cc, rv, nobs, max_iter, gtol, fopt, p_nagx_s_fs, &
                                            niter_tmp, conv_tmp)
                end if
                call nagarchx_variance_path(ret_cc, rv, p_nagx_s_fs%garch, h_nagx_fs(:,id))
                result%logL_full(id) = skewt_loglik(ret_cc(1:nobs), h_nagx_fs(1:nobs,id), &
                                                     p_nagx_s_fs%nu, p_nagx_s_fs%xi)

            case default
                result%logL_full(id) = result%nag_logL_full

            end select
            result%aic_full(id) = 2.0_dp*k - 2.0_dp*result%logL_full(id)
            result%bic_full(id) = logn_n*k  - 2.0_dp*result%logL_full(id)
        end do

        call print_comparison(result, dist_names)

        print '(/,A,F7.3,A,F7.3,A,F7.3)', &
              "read_sec=", read_sec, "  fit_sec=", fit_sec, "  elapsed_sec=", elapsed_sec

        deallocate(ret_cc, rv, h_nag, h_nagx, h_nag_fs, h_nagx_fs)
    end subroutine compare_one_file

    subroutine print_comparison(r, dist_names)
        type(nagarchx_dist_result_t), intent(in) :: r
        character(len=*),             intent(in) :: dist_names(:)
        integer :: ndist, id, best_id
        real(dp) :: best

        ndist = size(dist_names)

        ! Train/test table
        print '(/,A)', "Daily CC return forecast comparison (AIC/BIC from training logL)"
        print '(A,I0,A,I0)', "Training obs: ", r%ntrain, "  test obs: ", r%ntest
        print '(A)', repeat("-", 108)
        print '(A12,1X,A3,2(3X,A12),2(3X,A11),3X,A9,3X,A8,3X,A4,1X,A4)', &
              "Model", "k", "logL_train", "logL_test", "AIC", "BIC", "QLIKE", "vol_ann%", "iter", "conv"
        print '(A)', repeat("-", 108)
        print '(A12,1X,I3,2(3X,F12.3),2(3X,F11.3),3X,F9.5,3X,F8.3,3X,I4,1X,L1)', &
              "NAGARCH", k_nag, r%nag_logL_tr, r%nag_logL_te, r%nag_aic, r%nag_bic, &
              r%nag_qlike, r%nag_vol, r%nag_niter, r%nag_conv
        do id = 1, ndist
            print '(A12,1X,I3,2(3X,F12.3),2(3X,F11.3),3X,F9.5,3X,F8.3,3X,I4,1X,L1)', &
                  trim(model_label(dist_names(id))), dist_nparams(dist_names(id)), &
                  r%logL_tr(id), r%logL_te(id), r%aic(id), r%bic(id), &
                  r%qlike(id), r%vol(id), r%niter(id), r%conv(id)
        end do
        print '(A)', repeat("-", 108)

        best = r%nag_logL_te;  best_id = 0
        do id = 1, ndist
            if (r%logL_te(id) > best) then;  best = r%logL_te(id);  best_id = id;  end if
        end do
        if (best_id == 0) then
            print '(A)', "Best test logL: NAGARCH"
        else
            print '(A)', "Best test logL: "//trim(model_label(dist_names(best_id)))
        end if

        ! RV table
        print '(/,A)', "RV prediction (test set, h_t vs RV_t)"
        print '(A)', repeat("-", 46)
        print '(A12,3X,A9,3X,A10)', "Model", "QLIKE_RV", "Corr(h,RV)"
        print '(A)', repeat("-", 46)
        print '(A12,3X,F9.5,3X,F10.6)', "NAGARCH", r%nag_qlike_rv, r%nag_corr_rv
        do id = 1, ndist
            print '(A12,3X,F9.5,3X,F10.6)', trim(model_label(dist_names(id))), &
                  r%qlike_rv(id), r%corr_rv(id)
        end do
        print '(A)', repeat("-", 46)

        best = r%nag_qlike_rv;  best_id = 0
        do id = 1, ndist
            if (r%qlike_rv(id) < best) then;  best = r%qlike_rv(id);  best_id = id;  end if
        end do
        if (best_id == 0) then
            print '(A)', "Best QLIKE_RV: NAGARCH"
        else
            print '(A)', "Best QLIKE_RV: "//trim(model_label(dist_names(best_id)))
        end if

        ! Full-sample IC table
        call print_fullsample_comparison(r, dist_names)
    end subroutine print_comparison

    subroutine print_fullsample_comparison(r, dist_names)
        type(nagarchx_dist_result_t), intent(in) :: r
        character(len=*),             intent(in) :: dist_names(:)
        integer :: ndist, id, best_aic_id, best_bic_id
        real(dp) :: best_aic, best_bic

        ndist = size(dist_names)

        ! width = A12 + 1X + I3 + 3*(3X + F12.3) = 12+1+3+3*15 = 61 -> use 61
        print '(/,A,I0,A)', "Full-sample IC comparison (N=", r%nobs, ")"
        print '(A)', repeat("-", 61)
        print '(A12,1X,A3,3(3X,A12))', "Model", "k", "logL", "AIC", "BIC"
        print '(A)', repeat("-", 61)
        print '(A12,1X,I3,3(3X,F12.3))', &
              "NAGARCH", k_nag, r%nag_logL_full, r%nag_aic_full, r%nag_bic_full
        do id = 1, ndist
            print '(A12,1X,I3,3(3X,F12.3))', &
                  trim(model_label(dist_names(id))), dist_nparams(dist_names(id)), &
                  r%logL_full(id), r%aic_full(id), r%bic_full(id)
        end do
        print '(A)', repeat("-", 61)

        best_aic = r%nag_aic_full;  best_aic_id = 0
        best_bic = r%nag_bic_full;  best_bic_id = 0
        do id = 1, ndist
            if (r%aic_full(id) < best_aic) then;  best_aic = r%aic_full(id);  best_aic_id = id;  end if
            if (r%bic_full(id) < best_bic) then;  best_bic = r%bic_full(id);  best_bic_id = id;  end if
        end do
        if (best_aic_id == 0) then
            print '(A)', "Best AIC: NAGARCH"
        else
            print '(A)', "Best AIC: "//trim(model_label(dist_names(best_aic_id)))
        end if
        if (best_bic_id == 0) then
            print '(A)', "Best BIC: NAGARCH"
        else
            print '(A)', "Best BIC: "//trim(model_label(dist_names(best_bic_id)))
        end if
    end subroutine print_fullsample_comparison

    subroutine print_asset_summary(results, dist_names)
        type(nagarchx_dist_result_t), intent(in) :: results(:)
        character(len=*),             intent(in) :: dist_names(:)
        integer :: i, id, ndist, n, width, best_id
        integer :: wL_nag, wQ_nag, wRV_nag, wC_nag
        integer :: wA_nag, wB_nag   ! full-sample AIC, BIC
        integer, allocatable :: wL_nagx(:), wQ_nagx(:), wRV_nagx(:), wC_nagx(:)
        integer, allocatable :: wA_nagx(:), wB_nagx(:)
        real(dp) :: best
        character(len=14) :: abbrev

        n     = size(results)
        ndist = size(dist_names)
        allocate(wL_nagx(ndist), wQ_nagx(ndist), wRV_nagx(ndist), wC_nagx(ndist))
        allocate(wA_nagx(ndist), wB_nagx(ndist))
        wL_nag = 0;  wL_nagx = 0
        wQ_nag = 0;  wQ_nagx = 0
        wRV_nag = 0; wRV_nagx = 0
        wC_nag = 0;  wC_nagx = 0
        wA_nag = 0;  wA_nagx = 0
        wB_nag = 0;  wB_nagx = 0

        do i = 1, n
            best = results(i)%nag_logL_te
            do id = 1, ndist; best = max(best, results(i)%logL_te(id)); end do
            if (results(i)%nag_logL_te >= best) wL_nag = wL_nag + 1
            do id = 1, ndist
                if (results(i)%logL_te(id) >= best) wL_nagx(id) = wL_nagx(id) + 1
            end do

            best = results(i)%nag_qlike
            do id = 1, ndist; best = min(best, results(i)%qlike(id)); end do
            if (results(i)%nag_qlike <= best) wQ_nag = wQ_nag + 1
            do id = 1, ndist
                if (results(i)%qlike(id) <= best) wQ_nagx(id) = wQ_nagx(id) + 1
            end do

            best = results(i)%nag_qlike_rv
            do id = 1, ndist; best = min(best, results(i)%qlike_rv(id)); end do
            if (results(i)%nag_qlike_rv <= best) wRV_nag = wRV_nag + 1
            do id = 1, ndist
                if (results(i)%qlike_rv(id) <= best) wRV_nagx(id) = wRV_nagx(id) + 1
            end do

            best = results(i)%nag_corr_rv
            do id = 1, ndist; best = max(best, results(i)%corr_rv(id)); end do
            if (results(i)%nag_corr_rv >= best) wC_nag = wC_nag + 1
            do id = 1, ndist
                if (results(i)%corr_rv(id) >= best) wC_nagx(id) = wC_nagx(id) + 1
            end do

            best = results(i)%nag_aic_full
            do id = 1, ndist; best = min(best, results(i)%aic_full(id)); end do
            if (results(i)%nag_aic_full <= best) wA_nag = wA_nag + 1
            do id = 1, ndist
                if (results(i)%aic_full(id) <= best) wA_nagx(id) = wA_nagx(id) + 1
            end do

            best = results(i)%nag_bic_full
            do id = 1, ndist; best = min(best, results(i)%bic_full(id)); end do
            if (results(i)%nag_bic_full <= best) wB_nag = wB_nag + 1
            do id = 1, ndist
                if (results(i)%bic_full(id) <= best) wB_nagx(id) = wB_nagx(id) + 1
            end do
        end do

        print '(/,/,A,I0,A)', "Summary across ", n, " assets"

        ! width = 8 + (1+ndist)*(2+12) + (1+ndist)*(2+10) = 8 + (1+ndist)*26
        width = 8 + (1+ndist)*26

        ! CC return metrics
        print '(A)', "CC return forecasting metrics (test set)"
        print '(A)', repeat("-", width)
        write(*, '(A8)', advance='no') "Asset"
        write(*, '(2X,A12)', advance='no') "logL_te_NAG"
        do id = 1, ndist
            abbrev = "logL_te_"//trim(dist_abbrev(dist_names(id)))
            write(*, '(2X,A12)', advance='no') trim(abbrev)
        end do
        write(*, '(2X,A10)', advance='no') "QLIKE_NAG"
        do id = 1, ndist
            abbrev = "QLIKE_"//trim(dist_abbrev(dist_names(id)))
            write(*, '(2X,A10)', advance='no') trim(abbrev)
        end do
        write(*, '()')
        print '(A)', repeat("-", width)
        do i = 1, n
            write(*, '(A8)', advance='no') results(i)%asset
            write(*, '(2X,F12.3)', advance='no') results(i)%nag_logL_te
            do id = 1, ndist; write(*, '(2X,F12.3)', advance='no') results(i)%logL_te(id); end do
            write(*, '(2X,F10.5)', advance='no') results(i)%nag_qlike
            do id = 1, ndist; write(*, '(2X,F10.5)', advance='no') results(i)%qlike(id); end do
            write(*, '()')
        end do
        print '(A)', repeat("-", width)
        write(*, '(A,I0,A,I0)', advance='no') "Best logL_test: NAGARCH=", wL_nag, "/", n
        do id = 1, ndist
            write(*, '(2X,A,A,I0,A,I0)', advance='no') &
                  trim(model_label(dist_names(id))), "=", wL_nagx(id), "/", n
        end do
        write(*, '()')
        write(*, '(A,I0,A,I0)', advance='no') "Best QLIKE:     NAGARCH=", wQ_nag, "/", n
        do id = 1, ndist
            write(*, '(2X,A,A,I0,A,I0)', advance='no') &
                  trim(model_label(dist_names(id))), "=", wQ_nagx(id), "/", n
        end do
        write(*, '()')

        ! RV prediction metrics
        print '(/,A)', "RV prediction metrics (test set, h_t vs RV_t)"
        print '(A)', repeat("-", width)
        write(*, '(A8)', advance='no') "Asset"
        write(*, '(2X,A12)', advance='no') "QLIKERV_NAG"
        do id = 1, ndist
            abbrev = "QLIKERV_"//trim(dist_abbrev(dist_names(id)))
            write(*, '(2X,A12)', advance='no') trim(abbrev)
        end do
        write(*, '(2X,A10)', advance='no') "CorrRV_NAG"
        do id = 1, ndist
            abbrev = "CorrRV_"//trim(dist_abbrev(dist_names(id)))
            write(*, '(2X,A10)', advance='no') trim(abbrev)
        end do
        write(*, '()')
        print '(A)', repeat("-", width)
        do i = 1, n
            write(*, '(A8)', advance='no') results(i)%asset
            write(*, '(2X,F12.5)', advance='no') results(i)%nag_qlike_rv
            do id = 1, ndist; write(*, '(2X,F12.5)', advance='no') results(i)%qlike_rv(id); end do
            write(*, '(2X,F10.6)', advance='no') results(i)%nag_corr_rv
            do id = 1, ndist; write(*, '(2X,F10.6)', advance='no') results(i)%corr_rv(id); end do
            write(*, '()')
        end do
        print '(A)', repeat("-", width)
        write(*, '(A,I0,A,I0)', advance='no') "Best QLIKE_RV:  NAGARCH=", wRV_nag, "/", n
        do id = 1, ndist
            write(*, '(2X,A,A,I0,A,I0)', advance='no') &
                  trim(model_label(dist_names(id))), "=", wRV_nagx(id), "/", n
        end do
        write(*, '()')
        write(*, '(A,I0,A,I0)', advance='no') "Best Corr_RV:   NAGARCH=", wC_nag, "/", n
        do id = 1, ndist
            write(*, '(2X,A,A,I0,A,I0)', advance='no') &
                  trim(model_label(dist_names(id))), "=", wC_nagx(id), "/", n
        end do
        write(*, '()')

        ! Full-sample AIC/BIC
        ! width = 8 + (1+ndist)*(2+12) + (1+ndist)*(2+10) same formula
        print '(/,A)', "Full-sample IC (lower is better)"
        print '(A)', repeat("-", width)
        write(*, '(A8)', advance='no') "Asset"
        write(*, '(2X,A12)', advance='no') "AIC_NAG"
        do id = 1, ndist
            abbrev = "AIC_"//trim(dist_abbrev(dist_names(id)))
            write(*, '(2X,A12)', advance='no') trim(abbrev)
        end do
        write(*, '(2X,A10)', advance='no') "BIC_NAG"
        do id = 1, ndist
            abbrev = "BIC_"//trim(dist_abbrev(dist_names(id)))
            write(*, '(2X,A10)', advance='no') trim(abbrev)
        end do
        write(*, '()')
        print '(A)', repeat("-", width)
        do i = 1, n
            write(*, '(A8)', advance='no') results(i)%asset
            write(*, '(2X,F12.3)', advance='no') results(i)%nag_aic_full
            do id = 1, ndist; write(*, '(2X,F12.3)', advance='no') results(i)%aic_full(id); end do
            write(*, '(2X,F10.3)', advance='no') results(i)%nag_bic_full
            do id = 1, ndist; write(*, '(2X,F10.3)', advance='no') results(i)%bic_full(id); end do
            write(*, '()')
        end do
        print '(A)', repeat("-", width)
        write(*, '(A,I0,A,I0)', advance='no') "Best AIC:       NAGARCH=", wA_nag, "/", n
        do id = 1, ndist
            write(*, '(2X,A,A,I0,A,I0)', advance='no') &
                  trim(model_label(dist_names(id))), "=", wA_nagx(id), "/", n
        end do
        write(*, '()')
        write(*, '(A,I0,A,I0)', advance='no') "Best BIC:       NAGARCH=", wB_nag, "/", n
        do id = 1, ndist
            write(*, '(2X,A,A,I0,A,I0)', advance='no') &
                  trim(model_label(dist_names(id))), "=", wB_nagx(id), "/", n
        end do
        write(*, '()')

        ! suppress unused variable warning
        best_id = 0
        deallocate(wL_nagx, wQ_nagx, wRV_nagx, wC_nagx, wA_nagx, wB_nagx)
    end subroutine print_asset_summary

    ! Log-likelihood under standardised t(nu) innovations.
    real(dp) function t_loglik(y, h, nu)
        real(dp), intent(in) :: y(:), h(:), nu
        real(dp) :: hh(size(y))
        hh = max(h, min_var)
        t_loglik = sum(-0.5_dp*log(hh) + logpdf_std(y/sqrt(hh), dist_t, nu))
    end function t_loglik

    ! Log-likelihood under FS skewed-t(nu, xi) innovations.
    real(dp) function skewt_loglik(y, h, nu, xi)
        real(dp), intent(in) :: y(:), h(:), nu, xi
        real(dp) :: hh, z
        integer  :: i
        skewt_loglik = 0.0_dp
        do i = 1, size(y)
            hh = max(h(i), min_var)
            z  = y(i) / sqrt(hh)
            skewt_loglik = skewt_loglik - 0.5_dp*log(hh) + log(max(pdf_fs_skewt(z, nu, xi), min_pdf))
        end do
    end function skewt_loglik

    ! Number of free parameters for each supported distribution.
    integer function dist_nparams(dist_name)
        character(len=*), intent(in) :: dist_name
        select case (trim(dist_name))
        case ("T");        dist_nparams = 6
        case ("FS_SKEWT"); dist_nparams = 7
        case default;      dist_nparams = 5   ! NORMAL and fallback
        end select
    end function dist_nparams

    ! Short row label for the comparison table.
    pure function model_label(dist_name) result(label)
        character(len=*), intent(in) :: dist_name
        character(len=12) :: label
        select case (trim(dist_name))
        case ("T");        label = "NAGARCH-X T"
        case ("FS_SKEWT"); label = "NAGARCH-X S"
        case default;      label = "NAGARCH-X N"
        end select
    end function model_label

    ! One-letter abbreviation used in summary table column headers.
    pure function dist_abbrev(dist_name) result(abbrev)
        character(len=*), intent(in) :: dist_name
        character(len=4) :: abbrev
        select case (trim(dist_name))
        case ("T");        abbrev = "T"
        case ("FS_SKEWT"); abbrev = "S"
        case default;      abbrev = "N"
        end select
    end function dist_abbrev

    ! Full display name for the legend.
    pure function dist_fullname(dist_name) result(name)
        character(len=*), intent(in) :: dist_name
        character(len=20) :: name
        select case (trim(dist_name))
        case ("T");        name = "Student-t"
        case ("FS_SKEWT"); name = "FS skewed-t"
        case default;      name = "Normal"
        end select
    end function dist_fullname

end module nagarchx_dist_compare_mod

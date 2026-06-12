! Compare direct joint GARCH/distribution fits with two-step warm-started fits.
!
! Usage:
!   xfit_gen_garch_dist_warm_returns.exe [prices_file] [output_csv]
!
! For each asset and GARCH model, this program first fits normal-noise GARCH,
! fits standardized iid distributions to the resulting residuals, and then uses
! those estimates as starts for joint GARCH/distribution likelihood maximization.
! It also runs the usual direct joint fit so timings and likelihoods can be
! compared against the two-step screening approximation.

module garch_dist_warm_returns_mod
    use date_mod, only: print_program_header
    use kind_mod, only: dp
    use csv_mod, only: read_price_csv, print_price_sample_info
    use stats_mod, only: mean, variance
    use distributions_mod, only: std_dist_id_from_name => dist_id_from_name, fit_dist_std
    use garch_fit_dist_mod, only: garch_dist_fit_result_t, fit_garch_dist_model, &
        garch_dist_variance_path, model_param_count, dist_param_count, dist_nu_value, &
        dist_xi_value, dist_alpha_value
    implicit none
    private

    public :: run_garch_dist_warm_returns

    character(len=*), parameter :: default_prices_file = "spy_efa_eem_tlt_lqd.csv"
    real(dp), parameter :: trading_days = 252.0_dp
    real(dp), parameter :: min_h = 1.0e-12_dp
    integer, parameter :: max_assets = 3
    integer, parameter :: max_iter = 100
    real(dp), parameter :: gtol = 1.0e-7_dp
    character(len=16), parameter :: models(*) = [character(len=16) :: &
        "SYMM_GARCH", "NAGARCH"]
    character(len=10), parameter :: dists(*) = [character(len=10) :: &
        "NORMAL", "T", "SECH", "GED", "LAPLACE", "LOGISTIC", "NIG", "FS_SKEWT"]

contains

    subroutine run_garch_dist_warm_returns()
        integer, allocatable :: dates(:)
        character(len=32), allocatable :: col_names(:)
        real(dp), allocatable :: prices(:,:), ret(:), h(:), z(:)
        type(garch_dist_fit_result_t) :: normal_fit, direct_fit, warm_fit
        character(len=256) :: prices_file, output_csv
        integer :: nprices, ncols, nobs, icol, imod, idist, csv_unit
        integer :: clock_rate, clock_start, clock_end, fit_start, fit_end
        integer :: dist_iter, nparam
        real(dp) :: ret_mean, vol_ann, h_log_sum, dist_logl, two_logl, two_aic, two_bic
        real(dp) :: direct_logl, direct_aic, direct_bic, warm_logl, warm_aic, warm_bic
        real(dp) :: normal_sec, dist_sec, direct_sec, warm_sec, elapsed_s
        real(dp) :: shape, shape2
        logical :: write_csv, dist_converged

        call print_program_header("xfit_gen_garch_dist_warm_returns.f90")
        prices_file = default_prices_file
        output_csv = ""
        if (command_argument_count() >= 1) call get_command_argument(1, prices_file)
        if (command_argument_count() >= 2) call get_command_argument(2, output_csv)
        write_csv = len_trim(output_csv) > 0

        call system_clock(clock_start, clock_rate)
        call read_price_csv(prices_file, dates, col_names, prices, max_col=max_assets)
        nprices = size(prices, 1)
        ncols = size(prices, 2)
        nobs = nprices - 1
        allocate(ret(nobs), h(nobs), z(nobs))

        call print_price_sample_info(trim(prices_file), dates, ncols)
        write(*,'(A,I0)') "Maximum fit iterations: ", max_iter
        write(*,'(A,I0,A,I0)') "GARCH models: ", size(models), "  distributions: ", size(dists)

        if (write_csv) then
            open(newunit=csv_unit, file=output_csv, status="replace", action="write")
            write(csv_unit,'(A)') "asset,model,dist,vol_ann_pct,two_logL,two_AIC,two_BIC,joint_logL,joint_AIC," // &
                "joint_BIC,warm_logL,warm_AIC,warm_BIC,nu,xi,dist_alpha,two_sec,joint_sec,warm_joint_sec," // &
                "warm_total_sec,joint_iter,joint_conv,warm_iter,warm_conv,dist_iter,dist_conv"
        end if

        write(*,'(/,A)') "Joint GARCH/distribution fits: direct vs two-step warm starts"
        write(*,'(A)') "Model            Dist       Asset     two_logL     joint_logL      warm_logL      two_AIC    joint_AIC     warm_AIC   two_sec joint_sec  warm_sec warm_total  j_it j_cv  w_it w_cv"
        write(*,'(A)') repeat("-", 174)

        do icol = 1, ncols
            ret = log(prices(2:nprices,icol) / prices(1:nprices-1,icol))
            ret_mean = mean(ret)
            ret = ret - ret_mean
            vol_ann = sqrt(max(variance(ret), 0.0_dp)*trading_days)*100.0_dp

            do imod = 1, size(models)
                call system_clock(fit_start)
                call fit_garch_dist_model(trim(models(imod)), "NORMAL", ret, max_iter, gtol, normal_fit)
                call garch_dist_variance_path(trim(models(imod)), ret, normal_fit%params, h)
                h = max(h, min_h)
                z = ret / sqrt(h)
                h_log_sum = sum(log(h))
                call system_clock(fit_end)
                normal_sec = real(fit_end - fit_start, dp) / real(clock_rate, dp)

                do idist = 1, size(dists)
                    call system_clock(fit_start)
                    call fit_standardized_dist(trim(dists(idist)), z, shape, shape2, dist_logl, &
                                               dist_iter, dist_converged)
                    call system_clock(fit_end)
                    dist_sec = real(fit_end - fit_start, dp) / real(clock_rate, dp)

                    two_logl = dist_logl - 0.5_dp*h_log_sum
                    nparam = model_param_count(trim(models(imod))) + dist_param_count(trim(dists(idist)))
                    two_aic = 2.0_dp*real(nparam, dp) - 2.0_dp*two_logl
                    two_bic = log(real(nobs, dp))*real(nparam, dp) - 2.0_dp*two_logl

                    call system_clock(fit_start)
                    call fit_garch_dist_model(trim(models(imod)), trim(dists(idist)), ret, max_iter, gtol, direct_fit)
                    call system_clock(fit_end)
                    direct_sec = real(fit_end - fit_start, dp) / real(clock_rate, dp)

                    call system_clock(fit_start)
                    call fit_garch_dist_model(trim(models(imod)), trim(dists(idist)), ret, max_iter, gtol, warm_fit, &
                                              start_params=normal_fit%params, start_shape=shape, start_shape2=shape2)
                    call system_clock(fit_end)
                    warm_sec = real(fit_end - fit_start, dp) / real(clock_rate, dp)

                    direct_logl = -real(nobs, dp)*direct_fit%nll
                    direct_aic = 2.0_dp*real(nparam, dp) - 2.0_dp*direct_logl
                    direct_bic = log(real(nobs, dp))*real(nparam, dp) - 2.0_dp*direct_logl
                    warm_logl = -real(nobs, dp)*warm_fit%nll
                    warm_aic = 2.0_dp*real(nparam, dp) - 2.0_dp*warm_logl
                    warm_bic = log(real(nobs, dp))*real(nparam, dp) - 2.0_dp*warm_logl

                    write(*,'(A16,1X,A10,1X,A9,3F14.2,3F13.2,4F9.3,1X,I5,1X,L4,1X,I5,1X,L4)') &
                        trim(models(imod)), trim(dists(idist)), trim(col_names(icol)), two_logl, direct_logl, &
                        warm_logl, two_aic, direct_aic, warm_aic, normal_sec + dist_sec, direct_sec, &
                        warm_sec, normal_sec + dist_sec + warm_sec, direct_fit%niter, direct_fit%converged, &
                        warm_fit%niter, warm_fit%converged
                    if (write_csv) then
                        write(csv_unit,'(A,",",A,",",A,",",F12.6,9(",",F16.6),3(",",A),4(",",F12.6),3(",",I0,",",L1))') &
                            trim(col_names(icol)), trim(models(imod)), trim(dists(idist)), vol_ann, two_logl, &
                            two_aic, two_bic, direct_logl, direct_aic, direct_bic, warm_logl, warm_aic, warm_bic, &
                            dist_nu_value(trim(dists(idist)), shape), dist_xi_value(trim(dists(idist)), shape2), &
                            dist_alpha_value(trim(dists(idist)), shape), normal_sec + dist_sec, direct_sec, &
                            warm_sec, normal_sec + dist_sec + warm_sec, direct_fit%niter, direct_fit%converged, &
                            warm_fit%niter, warm_fit%converged, dist_iter, dist_converged
                    end if
                end do
            end do
        end do
        write(*,'(A)') repeat("-", 174)
        if (write_csv) close(csv_unit)

        call system_clock(clock_end)
        elapsed_s = real(clock_end - clock_start, dp) / real(clock_rate, dp)
        write(*,'(/,A,F10.3,A)') "Elapsed wall time: ", elapsed_s, " seconds"
    end subroutine run_garch_dist_warm_returns

    subroutine fit_standardized_dist(dist_name, z, shape, shape2, logl, niter, converged)
        character(len=*), intent(in) :: dist_name
        real(dp), intent(in) :: z(:)
        real(dp), intent(out) :: shape, shape2, logl
        integer, intent(out) :: niter
        logical, intent(out) :: converged
        integer :: dist_id

        dist_id = standardized_dist_id(dist_name)
        call fit_dist_std(z, size(z), dist_id, shape, logl, converged, niter, max_iter, xi_out=shape2)
    end subroutine fit_standardized_dist

    integer function standardized_dist_id(dist_name)
        character(len=*), intent(in) :: dist_name

        if (trim(dist_name) == "NIG") then
            standardized_dist_id = std_dist_id_from_name("NIG_SYM")
        else
            standardized_dist_id = std_dist_id_from_name(dist_name)
        end if
        if (standardized_dist_id == 0) error stop "standardized_dist_id: unsupported distribution"
    end function standardized_dist_id

end module garch_dist_warm_returns_mod

program xfit_gen_garch_dist_warm_returns
    use garch_dist_warm_returns_mod, only: run_garch_dist_warm_returns
    implicit none
    call run_garch_dist_warm_returns()
end program xfit_gen_garch_dist_warm_returns

! Fit the volatility models common to Python arch and this Fortran code.
! Reads daily prices, computes demeaned log returns, and fits ARCH(1),
! symmetric GARCH(1,1), QGARCH(1,1), FIGARCH(1,d,1), GJR-GARCH(1,1),
! EGARCH(1,1), APARCH(1,1), HARCH(1,5,22), RiskMetrics 2006, and symmetric/asymmetric MIDAS
! Hyperbolic with Normal errors.

program xfit_arch_common
    use kind_mod,       only: dp
    use math_const_mod, only: log_sqrt_2pi
    use csv_mod,        only: read_price_csv, print_price_sample_info
    use stats_mod,      only: mean
    use garch_types_mod, only: garch_params_t
    use garch_fit_mod,  only: fit_symm_garch, fit_qgarch, fit_figarch, fit_gjr_signed, fit_egarch, fit_aparch, fit_harch, &
                              fit_riskmetrics2006, &
                              fit_midas_hyperbolic, fit_midas_hyperbolic_asym, &
                              garch_skew_kurt, qgarch_skew_kurt, figarch_skew_kurt, gjr_skew_kurt, egarch_skew_kurt, aparch_skew_kurt, harch_skew_kurt, &
                              riskmetrics2006_skew_kurt, riskmetrics2006_variance, figarch_variance, midas_hyperbolic_skew_kurt, &
                              midas_hyperbolic_asym_skew_kurt, &
                              symm_garch_persist, qgarch_persist, figarch_persist, gjr_persist, egarch_persist, aparch_persist, harch_persist, &
                              riskmetrics2006_persist, midas_hyperbolic_persist, midas_hyperbolic_asym_persist, &
                              aparch_mean_variance, qgarch_mean_variance
    use bfgs_mod,       only: bfgs_minimize
    implicit none

    type fit_row_t
        character(len=16) :: model = ""
        character(len=16) :: asset = ""
        type(garch_params_t) :: params
        real(dp) :: persist = 0.0_dp
        real(dp) :: vol_ann = 0.0_dp
        real(dp) :: logl = 0.0_dp
        real(dp) :: aic = 0.0_dp
        real(dp) :: bic = 0.0_dp
        real(dp) :: skew = 0.0_dp
        real(dp) :: ekurt = 0.0_dp
        real(dp) :: seconds = 0.0_dp
        integer :: nparam = 0
        integer :: niter = 0
        integer :: aic_rank = 0
        integer :: bic_rank = 0
        logical :: converged = .false.
    end type fit_row_t

    character(len=*), parameter :: prices_file = "spy_efa_eem_tlt_lqd.csv"
    character(len=*), parameter :: csv_file = "arch/fortran_arch_common_results.csv"
    character(len=16), parameter :: models(*) = [character(len=16) :: &
        "ARCH1", "SYMM_GARCH", "QGARCH", "FIGARCH", "GJR_GARCH", "EGARCH", "APARCH", "HARCH", "RM2006", "MIDASHYP", &
        "MIDASHYP_ASYM"]
    real(dp), parameter :: trading_days = 252.0_dp
    integer,  parameter :: max_iter = 1000
    real(dp), parameter :: gtol = 1.0e-6_dp

    integer, allocatable :: dates(:)
    character(len=32), allocatable :: col_names(:)
    real(dp), allocatable :: prices(:,:), ret(:), arch_y(:)
    type(fit_row_t), allocatable :: rows(:)
    integer :: nprices, ncols, nobs, icol, imod, irow
    integer :: c0, c1, cr, t_read0, t_read1, t_fit0, t_fit1, t_print0, t_print1, t_write0, t_write1
    real(dp) :: read_s, fit_s, print_s, write_s, elapsed_s

    call system_clock(c0, cr)

    call system_clock(t_read0)
    call read_price_csv(prices_file, dates, col_names, prices)
    nprices = size(prices, 1)
    ncols = size(prices, 2)
    nobs = nprices - 1
    allocate(ret(nobs), arch_y(nobs), rows(ncols*size(models)))
    call system_clock(t_read1)

    call system_clock(t_fit0)
    irow = 0
    do icol = 1, ncols
        ret = log(prices(2:nprices,icol) / prices(1:nprices-1,icol))
        ret = ret - mean(ret)
        arch_y = ret

        do imod = 1, size(models)
            irow = irow + 1
            call fit_model(models(imod), trim(col_names(icol)), ret, rows(irow))
        end do
    end do
    call system_clock(t_fit1)

    call rank_results(rows, irow)

    call system_clock(t_print0)
    call print_results(rows, irow, dates, ncols)
    call system_clock(t_print1)

    call system_clock(t_write0)
    call write_results_csv(rows, irow)
    call system_clock(t_write1)

    call system_clock(c1)
    read_s = real(t_read1 - t_read0, dp) / real(cr, dp)
    fit_s = real(t_fit1 - t_fit0, dp) / real(cr, dp)
    print_s = real(t_print1 - t_print0, dp) / real(cr, dp)
    write_s = real(t_write1 - t_write0, dp) / real(cr, dp)
    elapsed_s = real(c1 - c0, dp) / real(cr, dp)
    write(*,'(/,A)') "Timing:"
    write(*,'(A,F10.3,A)') "  read/prep:    ", read_s, " seconds"
    write(*,'(A,F10.3,A)') "  fit models:   ", fit_s, " seconds"
    write(*,'(A,F10.3,A)') "  print output: ", print_s, " seconds"
    write(*,'(A,F10.3,A)') "  write CSV:    ", write_s, " seconds"
    write(*,'(A,F10.3,A)') "  elapsed wall: ", elapsed_s, " seconds"

contains

    subroutine fit_model(model_name, asset, y, row)
        character(len=*), intent(in) :: model_name, asset
        real(dp), intent(in) :: y(:)
        type(fit_row_t), intent(out) :: row
        integer :: t0, t1, rate
        real(dp) :: fopt, h_unc
        real(dp), allocatable :: variance(:)

        call system_clock(t0, rate)
        row = fit_row_t()
        row%model = model_name
        row%asset = asset

        select case (trim(model_name))
        case ("ARCH1")
            call fit_arch1(y, fopt, row%params, row%niter, row%converged)
            row%persist = row%params%alpha
            h_unc = row%params%omega / max(1.0_dp - row%persist, 1.0e-8_dp)
            row%nparam = 2
            call arch1_skew_kurt(y, row%params, row%skew, row%ekurt)
        case ("SYMM_GARCH")
            call fit_symm_garch(y, max_iter, gtol, fopt, row%params, row%niter, row%converged)
            row%persist = symm_garch_persist(row%params)
            h_unc = row%params%omega / max(1.0_dp - row%persist, 1.0e-8_dp)
            row%nparam = 3
            call garch_skew_kurt(y, row%params, row%skew, row%ekurt)
        case ("QGARCH", "QUADRATIC_GARCH")
            call fit_qgarch(y, max_iter, gtol, fopt, row%params, row%niter, row%converged)
            row%persist = qgarch_persist(row%params)
            h_unc = qgarch_mean_variance(row%params)
            row%nparam = 4
            call qgarch_skew_kurt(y, row%params, row%skew, row%ekurt)
            row%model = "QGARCH"
        case ("FIGARCH")
            call fit_figarch(y, max_iter, gtol, fopt, row%params, row%niter, row%converged)
            row%persist = figarch_persist(row%params)
            allocate(variance(size(y)))
            call figarch_variance(y, row%params, variance)
            h_unc = sum(variance) / real(size(y), dp)
            deallocate(variance)
            row%nparam = 4
            call figarch_skew_kurt(y, row%params, row%skew, row%ekurt)
        case ("GJR_GARCH")
            call fit_gjr_signed(y, max_iter, gtol, fopt, row%params, row%niter, row%converged)
            row%persist = gjr_persist(row%params)
            h_unc = row%params%omega / max(1.0_dp - row%persist, 1.0e-8_dp)
            row%nparam = 4
            call gjr_skew_kurt(y, row%params, row%skew, row%ekurt)
        case ("EGARCH")
            call fit_egarch(y, max_iter, gtol, fopt, row%params, row%niter, row%converged)
            row%persist = egarch_persist(row%params)
            h_unc = exp(row%params%omega / max(1.0_dp - row%params%beta, 1.0e-8_dp))
            row%nparam = 4
            call egarch_skew_kurt(y, row%params, row%skew, row%ekurt)
        case ("APARCH")
            call fit_aparch(y, max_iter, gtol, fopt, row%params, row%niter, row%converged)
            row%persist = aparch_persist(row%params)
            h_unc = aparch_mean_variance(y, row%params)
            row%nparam = 5
            call aparch_skew_kurt(y, row%params, row%skew, row%ekurt)
        case ("HARCH")
            call fit_harch(y, max_iter, gtol, fopt, row%params, row%niter, row%converged)
            row%persist = harch_persist(row%params)
            h_unc = row%params%omega / max(1.0_dp - row%persist, 1.0e-8_dp)
            row%nparam = 4
            call harch_skew_kurt(y, row%params, row%skew, row%ekurt)
        case ("RM2006", "RISKMETRICS2006")
            call fit_riskmetrics2006(y, fopt, row%params, row%niter, row%converged)
            row%persist = riskmetrics2006_persist()
            allocate(variance(size(y)))
            call riskmetrics2006_variance(y, variance)
            h_unc = sum(variance) / real(size(y), dp)
            deallocate(variance)
            row%nparam = 0
            call riskmetrics2006_skew_kurt(y, row%skew, row%ekurt)
            row%model = "RM2006"
        case ("MIDASHYP", "MIDAS_HYPERBOLIC")
            call fit_midas_hyperbolic(y, max_iter, gtol, fopt, row%params, row%niter, row%converged)
            row%persist = midas_hyperbolic_persist(row%params)
            h_unc = row%params%omega / max(1.0_dp - row%persist, 1.0e-8_dp)
            row%nparam = 3
            call midas_hyperbolic_skew_kurt(y, row%params, row%skew, row%ekurt)
            row%model = "MIDASHYP"
        case ("MIDASHYP_ASYM", "MIDAS_HYPERBOLIC_ASYM")
            call fit_midas_hyperbolic_asym(y, max_iter, gtol, fopt, row%params, row%niter, row%converged)
            row%persist = midas_hyperbolic_asym_persist(row%params)
            h_unc = row%params%omega / max(1.0_dp - row%persist, 1.0e-8_dp)
            row%nparam = 4
            call midas_hyperbolic_asym_skew_kurt(y, row%params, row%skew, row%ekurt)
            row%model = "MIDASHYP_ASYM"
        case default
            print '(A,A)', "Unknown model: ", trim(model_name)
            error stop
        end select

        row%vol_ann = sqrt(trading_days * h_unc) * 100.0_dp
        row%logl = -real(size(y), dp) * fopt
        row%aic = 2.0_dp * real(row%nparam, dp) - 2.0_dp * row%logl
        row%bic = log(real(size(y), dp)) * real(row%nparam, dp) - 2.0_dp * row%logl
        call system_clock(t1)
        row%seconds = real(t1 - t0, dp) / real(rate, dp)
    end subroutine fit_model

    subroutine rank_results(rows, nrow)
        type(fit_row_t), intent(inout) :: rows(:)
        integer, intent(in) :: nrow
        integer :: i, j

        do i = 1, nrow
            rows(i)%aic_rank = 1
            rows(i)%bic_rank = 1
            do j = 1, nrow
                if (trim(rows(j)%asset) /= trim(rows(i)%asset)) cycle
                if (rows(j)%aic < rows(i)%aic) rows(i)%aic_rank = rows(i)%aic_rank + 1
                if (rows(j)%bic < rows(i)%bic) rows(i)%bic_rank = rows(i)%bic_rank + 1
            end do
        end do
    end subroutine rank_results

    subroutine fit_arch1(y, f_best, params, niter_best, converged_best)
        real(dp), intent(in) :: y(:)
        real(dp), intent(out) :: f_best
        type(garch_params_t), intent(out) :: params
        integer, intent(out) :: niter_best
        logical, intent(out) :: converged_best
        integer, parameter :: np = 2, n_start = 4
        real(dp), parameter :: start_alpha(n_start) = [0.10_dp, 0.30_dp, 0.50_dp, 0.70_dp]
        real(dp) :: p(np), p0(np), p_best(np), f_try, omega0
        integer :: i, niter_try
        logical :: converged_try

        arch_y = y
        f_best = huge(1.0_dp)
        p_best = 0.0_dp
        niter_best = 0
        converged_best = .false.
        do i = 1, n_start
            omega0 = max((1.0_dp - start_alpha(i))*sum(y**2)/real(size(y), dp), 1.0e-12_dp)
            call arch1_inv_transform(omega0, start_alpha(i), p0)
            p = p0
            call bfgs_minimize(arch1_obj, p, np, max_iter, gtol, f_try, niter_try, converged_try)
            if (f_try < f_best) then
                f_best = f_try
                p_best = p
                niter_best = niter_try
                converged_best = converged_try
            end if
        end do
        params = garch_params_t()
        call arch1_transform(p_best, params%omega, params%alpha)
    end subroutine fit_arch1

    subroutine arch1_obj(p, np, f, g)
        integer, intent(in) :: np
        real(dp), intent(in) :: p(np)
        real(dp), intent(out) :: f, g(np)
        real(dp) :: pp(np), pm(np), fp, fm, step
        integer :: i

        f = arch1_value(p)
        do i = 1, np
            step = 1.0e-5_dp * max(1.0_dp, abs(p(i)))
            pp = p
            pm = p
            pp(i) = pp(i) + step
            pm(i) = pm(i) - step
            fp = arch1_value(pp)
            fm = arch1_value(pm)
            g(i) = (fp - fm) / (2.0_dp*step)
        end do
    end subroutine arch1_obj

    real(dp) function arch1_value(p)
        real(dp), intent(in) :: p(:)
        real(dp) :: omega, alpha, h
        integer :: t

        call arch1_transform(p, omega, alpha)
        h = max(omega / max(1.0_dp - alpha, 1.0e-8_dp), sum(arch_y**2)/real(size(arch_y), dp), 1.0e-12_dp)
        arch1_value = 0.0_dp
        do t = 1, size(arch_y)
            h = max(h, 1.0e-12_dp)
            arch1_value = arch1_value + log_sqrt_2pi + 0.5_dp*(log(h) + arch_y(t)**2/h)
            h = omega + alpha*arch_y(t)**2
        end do
        arch1_value = arch1_value / real(size(arch_y), dp)
    end function arch1_value

    subroutine arch1_transform(p, omega, alpha)
        real(dp), intent(in) :: p(:)
        real(dp), intent(out) :: omega, alpha

        omega = exp(p(1))
        alpha = 1.0_dp / (1.0_dp + exp(-p(2)))
        alpha = min(max(alpha, 1.0e-8_dp), 1.0_dp - 1.0e-8_dp)
    end subroutine arch1_transform

    subroutine arch1_inv_transform(omega, alpha, p)
        real(dp), intent(in) :: omega, alpha
        real(dp), intent(out) :: p(:)
        real(dp) :: a

        a = min(max(alpha, 1.0e-8_dp), 1.0_dp - 1.0e-8_dp)
        p(1) = log(max(omega, 1.0e-12_dp))
        p(2) = log(a / (1.0_dp - a))
    end subroutine arch1_inv_transform

    subroutine arch1_skew_kurt(y, params, skew, ekurt)
        real(dp), intent(in) :: y(:)
        type(garch_params_t), intent(in) :: params
        real(dp), intent(out) :: skew, ekurt
        real(dp), allocatable :: z(:)
        real(dp) :: h
        integer :: t

        allocate(z(size(y)))
        h = max(params%omega / max(1.0_dp - params%alpha, 1.0e-8_dp), sum(y**2)/real(size(y), dp), 1.0e-12_dp)
        do t = 1, size(y)
            h = max(h, 1.0e-12_dp)
            z(t) = y(t) / sqrt(h)
            h = params%omega + params%alpha*y(t)**2
        end do
        call moments(z, skew, ekurt)
        deallocate(z)
    end subroutine arch1_skew_kurt

    subroutine moments(z, skew, ekurt)
        real(dp), intent(in) :: z(:)
        real(dp), intent(out) :: skew, ekurt
        real(dp) :: avg, m2, m3, m4, dz, rn
        integer :: i

        rn = real(size(z), dp)
        avg = sum(z) / rn
        m2 = 0.0_dp
        m3 = 0.0_dp
        m4 = 0.0_dp
        do i = 1, size(z)
            dz = z(i) - avg
            m2 = m2 + dz**2
            m3 = m3 + dz**3
            m4 = m4 + dz**4
        end do
        m2 = m2 / rn
        if (m2 <= 0.0_dp) then
            skew = 0.0_dp
            ekurt = 0.0_dp
        else
            skew = (m3 / rn) / m2**1.5_dp
            ekurt = (m4 / rn) / m2**2 - 3.0_dp
        end if
    end subroutine moments

    subroutine print_results(rows, nrow, dates, nassets)
        type(fit_row_t), intent(in) :: rows(:)
        integer, intent(in) :: nrow, dates(:), nassets
        integer :: i

        call print_price_sample_info(prices_file, dates, nassets)
        write(*,'(A)') "Package: Fortran GARCH-BFGS"
        write(*,'(A)') "Model            Asset        omega   alpha   gamma    beta   delta  persist  vol_ann%        logL         AIC         BIC #param iter conv    skew   ekurt      sec AIC_rank BIC_rank"
        write(*,'(A)') repeat("-", 192)
        do i = 1, nrow
            write(*,'(A16,1X,A9,ES12.3,5F8.4,F10.2,3F12.2,I7,I5,1X,L1,2F9.3,F9.3,2I9)') &
                trim(rows(i)%model), trim(rows(i)%asset), rows(i)%params%omega, rows(i)%params%alpha, &
                rows(i)%params%gamma, rows(i)%params%beta, rows(i)%params%theta, rows(i)%persist, rows(i)%vol_ann, &
                rows(i)%logl, rows(i)%aic, rows(i)%bic, rows(i)%nparam, rows(i)%niter, rows(i)%converged, &
                rows(i)%skew, rows(i)%ekurt, rows(i)%seconds, rows(i)%aic_rank, rows(i)%bic_rank
        end do
        write(*,'(/,A,A)') "Wrote ", trim(csv_file)
    end subroutine print_results

    subroutine write_results_csv(rows, nrow)
        type(fit_row_t), intent(in) :: rows(:)
        integer, intent(in) :: nrow
        integer :: unit, i

        open(newunit=unit, file=csv_file, status="replace", action="write")
        write(unit,'(A)') "asset,model,omega,alpha,gamma,beta,delta,persist,vol_ann_pct,logL,AIC,BIC,nparam,iter,conv,skew,ekurt,sec,AIC_rank,BIC_rank"
        do i = 1, nrow
            write(unit,'(A,",",A,",",ES24.16,",",ES24.16,",",ES24.16,",",ES24.16,",",ES24.16,",",ES24.16,",",ES24.16,",",ES24.16,",",ES24.16,",",ES24.16,",",I0,",",I0,",",L1,",",ES24.16,",",ES24.16,",",ES24.16,",",I0,",",I0)') &
                trim(rows(i)%asset), trim(rows(i)%model), rows(i)%params%omega, rows(i)%params%alpha, &
                rows(i)%params%gamma, rows(i)%params%beta, rows(i)%params%theta, rows(i)%persist, rows(i)%vol_ann, rows(i)%logl, &
                rows(i)%aic, rows(i)%bic, rows(i)%nparam, rows(i)%niter, rows(i)%converged, &
                rows(i)%skew, rows(i)%ekurt, rows(i)%seconds, rows(i)%aic_rank, rows(i)%bic_rank
        end do
        close(unit)
    end subroutine write_results_csv

end program xfit_arch_common

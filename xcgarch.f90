! Fit Component Symmetric GARCH (CSGARCH / Engle-Lee 1999) and optionally
! Component NAGARCH (CNAGARCH, single-theta asymmetric extension) and GARCH(1,1)
! to an ETF return panel.
!
! CSGARCH:
!   q_t = omega + rho*q_{t-1} + phi*(y_{t-1}^2 - h_{t-1})
!   h_t = q_t + alpha*(y_{t-1}^2 - q_{t-1}) + beta*(h_{t-1} - q_{t-1})
!
! CNAGARCH (phi=0 nests NAGARCH; theta=0 nests CSGARCH; both nest GARCH):
!   q_t = omega + rho*q_{t-1} + phi*((y_{t-1} - theta*sqrt(h_{t-1}))^2 - h_{t-1})
!   h_t = q_t + alpha*((y_{t-1} - theta*sqrt(h_{t-1}))^2 - q_{t-1}) + beta*(h_{t-1} - q_{t-1})
!
! Usage: xcgarch.exe [prices_file] [--no-csv]

program xcgarch
    use date_mod,        only: print_program_header, date_label
    use kind_mod,        only: dp
    use math_const_mod,  only: log_sqrt_2pi
    use csv_mod,         only: read_price_csv, print_price_sample_info
    use stats_mod,       only: mean, column_summary_stats, print_column_summary_stats
    use garch_types_mod, only: garch_params_t
    use garch_fit_mod,   only: fit_csgarch,    csgarch_persist,    csgarch_skew_kurt, &
                                fit_symm_garch, symm_garch_persist, garch_skew_kurt
    use bfgs_mod,        only: bfgs_minimize
    implicit none

    logical, parameter :: print_ret_stats = .true.
    logical, parameter :: fit_garch_1_1  = .true.
    logical, parameter :: fit_cnagarch   = .true.

    character(len=*), parameter :: default_prices_file = "spy_efa_eem_tlt_lqd.csv"
    integer,  parameter :: max_iter     = 500
    integer,  parameter :: csgarch_np   = 5
    integer,  parameter :: garch_np     = 3
    integer,  parameter :: cnagarch_np  = 6
    integer,  parameter :: trading_days = 252
    real(dp), parameter :: gtol         = 1.0e-7_dp

    integer,           allocatable :: dates(:)
    character(len=32), allocatable :: col_names(:)
    real(dp),          allocatable :: prices(:,:), ret(:), h(:), q(:), h_g(:)
    real(dp),          allocatable :: cnagarch_data(:)   ! host-associated by objective

    type(garch_params_t), allocatable :: all_cp(:), all_gp(:), all_np(:)
    real(dp), allocatable :: all_cf(:), all_gf(:), all_nf(:)
    integer,  allocatable :: all_cn(:), all_gn(:), all_nn(:)
    logical,  allocatable :: all_cc(:), all_gc(:), all_nc(:)

    integer  :: nprices, ncols, nobs, icol, t, iarg, csv_unit
    real(dp) :: logl, aic, bic, logl_g, aic_g, bic_g, logl_n, aic_n, bic_n
    real(dp) :: skew, ekurt, persist_q, vol_ann
    integer  :: clock_start, clock_end, clock_rate
    real(dp) :: elapsed_s
    logical  :: write_csv
    character(len=256) :: prices_file, arg, csv_file

    call print_program_header("xcgarch.f90")

    prices_file = default_prices_file
    write_csv   = .true.
    do iarg = 1, command_argument_count()
        call get_command_argument(iarg, arg)
        if (trim(arg) == "--no-csv") then
            write_csv = .false.
        else
            prices_file = arg
        end if
    end do

    call system_clock(clock_start, clock_rate)
    call read_price_csv(prices_file, dates, col_names, prices)
    nprices = size(prices, 1)
    ncols   = size(prices, 2)
    nobs    = nprices - 1
    allocate(ret(nobs), h(nobs), q(nobs), h_g(nobs), cnagarch_data(nobs))
    allocate(all_cp(ncols), all_cf(ncols), all_cn(ncols), all_cc(ncols))
    if (fit_garch_1_1) &
        allocate(all_gp(ncols), all_gf(ncols), all_gn(ncols), all_gc(ncols))
    if (fit_cnagarch) &
        allocate(all_np(ncols), all_nf(ncols), all_nn(ncols), all_nc(ncols))

    call print_price_sample_info(prices_file, dates, ncols)

    write(*,'(A)') "CSGARCH:  q_t = omega + rho*q_{t-1} + phi*(y^2_{t-1} - h_{t-1})"
    write(*,'(A)') "          h_t = q_t + alpha*(y^2_{t-1}-q_{t-1}) + beta*(h_{t-1}-q_{t-1})"
    if (fit_cnagarch) then
        write(*,'(A)') "CNAGARCH: replace y^2_{t-1} with (y_{t-1} - theta*sqrt(h_{t-1}))^2"
    end if

    ! --- return summary statistics ---
    if (print_ret_stats) then
        block
            real(dp), allocatable :: all_ret(:,:), med(:), avg(:), sdev(:), sk(:), ek(:), xmn(:), xmx(:)
            integer :: j
            allocate(all_ret(nobs, ncols))
            do j = 1, ncols
                all_ret(:,j) = log(prices(2:nprices,j) / prices(1:nprices-1,j)) * 100.0_dp
            end do
            call column_summary_stats(all_ret, med, avg, sdev, sk, ek, xmn, xmx)
            avg  = avg  * real(trading_days, dp)
            sdev = sdev * sqrt(real(trading_days, dp))
            call print_column_summary_stats( &
                "Log-return statistics (mean and sd annualized, other stats daily %)", &
                col_names, med, avg, sdev, sk, ek, xmn, xmx)
        end block
    end if

    ! --- fit all assets ---
    do icol = 1, ncols
        ret            = log(prices(2:nprices,icol) / prices(1:nprices-1,icol))
        ret            = ret - mean(ret)
        cnagarch_data  = ret
        call fit_csgarch(ret, max_iter, gtol, all_cf(icol), all_cp(icol), &
                         all_cn(icol), all_cc(icol))
        if (fit_garch_1_1) &
            call fit_symm_garch(ret, max_iter, gtol, all_gf(icol), all_gp(icol), &
                                all_gn(icol), all_gc(icol))
        if (fit_cnagarch) &
            call fit_cnagarch_local(all_cp(icol), all_nf(icol), all_np(icol), &
                                    all_nn(icol), all_nc(icol))
    end do

    ! --- CSGARCH parameter table ---
    write(*,'(/,A)') repeat("-", 134)
    write(*,'(A)') &
        "CSGARCH   " // &
        "     omega   alpha    beta     rho     phi  persist  unc_vol%" // &
        "       logL        AIC        BIC  iter conv    skew   ekurt"
    write(*,'(A)') repeat("-", 134)
    do icol = 1, ncols
        ret      = log(prices(2:nprices,icol) / prices(1:nprices-1,icol))
        ret      = ret - mean(ret)
        call csgarch_hq(ret, all_cp(icol), h, q)
        call csgarch_skew_kurt(ret, all_cp(icol), skew, ekurt)
        persist_q = csgarch_persist(all_cp(icol))
        vol_ann   = sqrt(real(trading_days, dp) * sum(h) / real(nobs, dp)) * 100.0_dp
        logl      = -real(nobs, dp) * all_cf(icol)
        aic       = 2.0_dp * real(csgarch_np, dp) - 2.0_dp * logl
        bic       = log(real(nobs, dp)) * real(csgarch_np, dp) - 2.0_dp * logl
        write(*,'(A10,ES12.3,5F8.4,F9.2,3F11.2,I6,1X,L1,2F8.3)') &
            trim(col_names(icol)), all_cp(icol)%omega, &
            all_cp(icol)%alpha, all_cp(icol)%beta, &
            all_cp(icol)%extra1, all_cp(icol)%extra2, &
            persist_q, vol_ann, logl, aic, bic, &
            all_cn(icol), all_cc(icol), skew, ekurt
    end do
    write(*,'(A)') repeat("-", 134)

    ! --- CNAGARCH parameter table ---
    if (fit_cnagarch) then
        write(*,'(/,A)') repeat("-", 142)
        write(*,'(A)') &
            "CNAGARCH  " // &
            "     omega   alpha    beta     rho     phi   theta  t_pers  unc_vol%" // &
            "       logL        AIC        BIC  iter conv    skew   ekurt"
        write(*,'(A)') repeat("-", 142)
        do icol = 1, ncols
            ret      = log(prices(2:nprices,icol) / prices(1:nprices-1,icol))
            ret      = ret - mean(ret)
            call cnagarch_hq(ret, all_np(icol), h, q)
            call cnagarch_skew_kurt(ret, all_np(icol), skew, ekurt)
            vol_ann  = sqrt(real(trading_days, dp) * sum(h) / real(nobs, dp)) * 100.0_dp
            logl_n   = -real(nobs, dp) * all_nf(icol)
            aic_n    = 2.0_dp * real(cnagarch_np, dp) - 2.0_dp * logl_n
            bic_n    = log(real(nobs, dp)) * real(cnagarch_np, dp) - 2.0_dp * logl_n
            write(*,'(A10,ES12.3,6F8.4,F9.2,3F11.2,I6,1X,L1,2F8.3)') &
                trim(col_names(icol)), all_np(icol)%omega, &
                all_np(icol)%alpha, all_np(icol)%beta, &
                all_np(icol)%extra1, all_np(icol)%extra2, all_np(icol)%theta, &
                all_np(icol)%alpha*(1.0_dp + all_np(icol)%theta**2) + all_np(icol)%beta, &
                vol_ann, logl_n, aic_n, bic_n, &
                all_nn(icol), all_nc(icol), skew, ekurt
        end do
        write(*,'(A)') repeat("-", 142)
    end if

    ! --- GARCH(1,1) parameter table ---
    if (fit_garch_1_1) then
        write(*,'(/,A)') repeat("-", 114)
        write(*,'(A)') &
            "GARCH(1,1)" // &
            "       omega   alpha    beta  persist  unc_vol%" // &
            "       logL        AIC        BIC  iter conv    skew   ekurt"
        write(*,'(A)') repeat("-", 114)
        do icol = 1, ncols
            ret      = log(prices(2:nprices,icol) / prices(1:nprices-1,icol))
            ret      = ret - mean(ret)
            call garch11_h(ret, all_gp(icol), h_g)
            call garch_skew_kurt(ret, all_gp(icol), skew, ekurt)
            persist_q = symm_garch_persist(all_gp(icol))
            vol_ann   = sqrt(real(trading_days, dp) * sum(h_g) / real(nobs, dp)) * 100.0_dp
            logl_g    = -real(nobs, dp) * all_gf(icol)
            aic_g     = 2.0_dp * real(garch_np, dp) - 2.0_dp * logl_g
            bic_g     = log(real(nobs, dp)) * real(garch_np, dp) - 2.0_dp * logl_g
            write(*,'(A10,ES12.3,3F8.4,F9.2,3F11.2,I6,1X,L1,2F8.3)') &
                trim(col_names(icol)), all_gp(icol)%omega, &
                all_gp(icol)%alpha, all_gp(icol)%beta, &
                persist_q, vol_ann, logl_g, aic_g, bic_g, &
                all_gn(icol), all_gc(icol), skew, ekurt
        end do
        write(*,'(A)') repeat("-", 114)
    end if

    ! --- IC comparison table ---
    if (fit_garch_1_1 .or. fit_cnagarch) call print_ic_comparison()

    ! --- component volatility time-series: one CSV file per asset (CSGARCH) ---
    if (write_csv) then
        do icol = 1, ncols
            ret = log(prices(2:nprices,icol) / prices(1:nprices-1,icol))
            ret = ret - mean(ret)
            call csgarch_hq(ret, all_cp(icol), h, q)
            csv_file = "cgarch_" // trim(col_names(icol)) // ".csv"
            open(newunit=csv_unit, file=trim(csv_file), status="replace", action="write")
            write(csv_unit,'(A)') "date,ret_pct,tot_vol_pct,perm_vol_pct,tran_vol_pct,tran_frac"
            do t = 1, nobs
                write(csv_unit,'(A,5(",",F12.8))') &
                    date_label(dates(t+1)), &
                    ret(t) * 100.0_dp, &
                    sqrt(real(trading_days,dp) * h(t)) * 100.0_dp, &
                    sqrt(real(trading_days,dp) * q(t)) * 100.0_dp, &
                    sqrt(real(trading_days,dp) * max(h(t) - q(t), 0.0_dp)) * 100.0_dp, &
                    (h(t) - q(t)) / max(h(t), 1.0e-12_dp)
            end do
            close(csv_unit)
            write(*,'(/,A,A)') "Written: ", trim(csv_file)
        end do
    end if

    call system_clock(clock_end)
    elapsed_s = real(clock_end - clock_start, dp) / real(clock_rate, dp)
    write(*,'(/,A,F10.3,A)') "Elapsed wall time:", elapsed_s, " seconds"

    deallocate(dates, col_names, prices, ret, h, q, h_g, cnagarch_data)
    deallocate(all_cp, all_cf, all_cn, all_cc)
    if (fit_garch_1_1) deallocate(all_gp, all_gf, all_gn, all_gc)
    if (fit_cnagarch)  deallocate(all_np, all_nf, all_nn, all_nc)

contains

    ! -------------------------------------------------------------------------
    ! IC comparison table: shows AIC and BIC for all active models per asset,
    ! averages, and win counts.
    ! -------------------------------------------------------------------------
    subroutine print_ic_comparison()
        real(dp) :: a_c, b_c, a_n, b_n, a_g, b_g
        real(dp) :: sum_da_nc, sum_db_nc, sum_da_gc, sum_db_gc
        integer  :: win_c_a, win_n_a, win_g_a, win_c_b, win_n_b, win_g_b
        integer  :: ic, best_a, best_b
        real(dp) :: best_aic, best_bic

        sum_da_nc = 0.0_dp;  sum_db_nc = 0.0_dp
        sum_da_gc = 0.0_dp;  sum_db_gc = 0.0_dp
        win_c_a = 0;  win_n_a = 0;  win_g_a = 0
        win_c_b = 0;  win_n_b = 0;  win_g_b = 0

        write(*,'(/,A)') "IC comparison (lower is better)"
        write(*,'(A)') repeat("-", 98)
        write(*,'(A10)', advance='no') "Asset"
        write(*,'(A11)', advance='no') "AIC_CGARCH"
        if (fit_cnagarch)  write(*,'(A11)', advance='no') "AIC_CNAG"
        if (fit_garch_1_1) write(*,'(A11)', advance='no') "AIC_GARCH"
        write(*,'(A11)', advance='no') "BIC_CGARCH"
        if (fit_cnagarch)  write(*,'(A11)', advance='no') "BIC_CNAG"
        if (fit_garch_1_1) write(*,'(A11)', advance='no') "BIC_GARCH"
        write(*,*)
        write(*,'(A)') repeat("-", 98)

        do ic = 1, ncols
            a_c = 2.0_dp*real(csgarch_np,dp)  - 2.0_dp*(-real(nobs,dp)*all_cf(ic))
            b_c = log(real(nobs,dp))*real(csgarch_np,dp) - 2.0_dp*(-real(nobs,dp)*all_cf(ic))
            a_n = 0.0_dp;  b_n = 0.0_dp
            a_g = 0.0_dp;  b_g = 0.0_dp
            if (fit_cnagarch) then
                a_n = 2.0_dp*real(cnagarch_np,dp)  - 2.0_dp*(-real(nobs,dp)*all_nf(ic))
                b_n = log(real(nobs,dp))*real(cnagarch_np,dp) - 2.0_dp*(-real(nobs,dp)*all_nf(ic))
                sum_da_nc = sum_da_nc + (a_n - a_c)
                sum_db_nc = sum_db_nc + (b_n - b_c)
            end if
            if (fit_garch_1_1) then
                a_g = 2.0_dp*real(garch_np,dp)  - 2.0_dp*(-real(nobs,dp)*all_gf(ic))
                b_g = log(real(nobs,dp))*real(garch_np,dp) - 2.0_dp*(-real(nobs,dp)*all_gf(ic))
                sum_da_gc = sum_da_gc + (a_g - a_c)
                sum_db_gc = sum_db_gc + (b_g - b_c)
            end if

            ! find best model by AIC
            best_aic = a_c;  best_a = 1
            if (fit_cnagarch  .and. a_n < best_aic) then;  best_aic = a_n;  best_a = 2;  end if
            if (fit_garch_1_1 .and. a_g < best_aic) then;                   best_a = 3;  end if
            if (best_a == 1) win_c_a = win_c_a + 1
            if (best_a == 2) win_n_a = win_n_a + 1
            if (best_a == 3) win_g_a = win_g_a + 1

            ! find best model by BIC
            best_bic = b_c;  best_b = 1
            if (fit_cnagarch  .and. b_n < best_bic) then;  best_bic = b_n;  best_b = 2;  end if
            if (fit_garch_1_1 .and. b_g < best_bic) then;                   best_b = 3;  end if
            if (best_b == 1) win_c_b = win_c_b + 1
            if (best_b == 2) win_n_b = win_n_b + 1
            if (best_b == 3) win_g_b = win_g_b + 1

            write(*,'(A10,F11.2)', advance='no') trim(col_names(ic)), a_c
            if (fit_cnagarch)  write(*,'(F11.2)', advance='no') a_n
            if (fit_garch_1_1) write(*,'(F11.2)', advance='no') a_g
            write(*,'(F11.2)', advance='no') b_c
            if (fit_cnagarch)  write(*,'(F11.2)', advance='no') b_n
            if (fit_garch_1_1) write(*,'(F11.2)', advance='no') b_g
            write(*,*)
        end do

        write(*,'(A)') repeat("-", 98)

        ! averages row
        write(*,'(A10,A11)', advance='no') "Average", ""
        if (fit_cnagarch)  write(*,'(F11.2)', advance='no') sum_da_nc / ncols
        if (fit_garch_1_1) write(*,'(F11.2)', advance='no') sum_da_gc / ncols
        write(*,'(A11)',   advance='no') ""
        if (fit_cnagarch)  write(*,'(F11.2)', advance='no') sum_db_nc / ncols
        if (fit_garch_1_1) write(*,'(F11.2)', advance='no') sum_db_gc / ncols
        write(*,*)
        write(*,'(A)') "  (averages are differences from CGARCH; negative favors the other model)"

        write(*,'(A)') repeat("-", 98)

        ! win counts
        write(*,'(A,I0,A,I0,A,I0)') &
            "AIC best: CGARCH=", win_c_a, "  CNAGARCH=", win_n_a, "  GARCH=", win_g_a
        write(*,'(A,I0,A,I0,A,I0)') &
            "BIC best: CGARCH=", win_c_b, "  CNAGARCH=", win_n_b, "  GARCH=", win_g_b
    end subroutine print_ic_comparison

    ! -------------------------------------------------------------------------
    ! CSGARCH h/q decomposition
    ! -------------------------------------------------------------------------
    subroutine csgarch_hq(y, params, h_out, q_out)
        real(dp),             intent(in)  :: y(:)
        type(garch_params_t), intent(in)  :: params
        real(dp),             intent(out) :: h_out(:), q_out(:)
        real(dp) :: backcast, h_prev, q_prev, h, qq, shock
        integer  :: t

        backcast = max(sum(y**2) / real(size(y), dp), 1.0e-12_dp)
        q_prev   = params%omega / max(1.0_dp - params%extra1, 1.0e-8_dp)
        q_prev   = max(q_prev, 1.0e-12_dp)
        h_prev   = max(q_prev, backcast)
        shock    = backcast
        do t = 1, size(y)
            qq       = params%omega + params%extra1*q_prev + params%extra2*(shock - h_prev)
            qq       = max(qq, 1.0e-12_dp)
            h        = qq + params%alpha*(shock - q_prev) + params%beta*(h_prev - q_prev)
            h        = max(h, 1.0e-12_dp)
            h_out(t) = h;  q_out(t) = qq
            q_prev = qq;   h_prev = h;   shock = y(t)**2
        end do
    end subroutine csgarch_hq

    ! -------------------------------------------------------------------------
    ! GARCH(1,1) variance path
    ! -------------------------------------------------------------------------
    subroutine garch11_h(y, params, h_out)
        real(dp),             intent(in)  :: y(:)
        type(garch_params_t), intent(in)  :: params
        real(dp),             intent(out) :: h_out(:)
        real(dp) :: h, shock
        integer  :: t

        h     = max(sum(y**2) / real(size(y), dp), 1.0e-12_dp)
        shock = h
        do t = 1, size(y)
            h        = params%omega + params%alpha * shock + params%beta * h
            h        = max(h, 1.0e-12_dp)
            h_out(t) = h;   shock = y(t)**2
        end do
    end subroutine garch11_h

    ! -------------------------------------------------------------------------
    ! CNAGARCH parameter transforms
    ! -------------------------------------------------------------------------
    subroutine cnagarch_transform(p, omega, alpha, beta, rho, phi, theta)
        real(dp), intent(in)  :: p(cnagarch_np)
        real(dp), intent(out) :: omega, alpha, beta, rho, phi, theta
        real(dp) :: ea, eb, denom
        omega = exp(p(1))
        ea    = exp(p(2));  eb = exp(p(3));  denom = 1.0_dp + ea + eb
        alpha = ea / denom;  beta = eb / denom
        rho   = 1.0_dp   / (1.0_dp + exp(-p(4)))
        phi   = 0.25_dp  / (1.0_dp + exp(-p(5)))
        theta = p(6)
    end subroutine cnagarch_transform

    subroutine cnagarch_inv_transform(omega, alpha, beta, rho, phi, theta, p)
        real(dp), intent(in)  :: omega, alpha, beta, rho, phi, theta
        real(dp), intent(out) :: p(cnagarch_np)
        real(dp) :: slack, rr, pp
        slack = max(1.0_dp - alpha - beta, 1.0e-8_dp)
        rr    = min(max(rho,       1.0e-8_dp), 1.0_dp - 1.0e-8_dp)
        pp    = min(max(phi/0.25_dp, 1.0e-8_dp), 1.0_dp - 1.0e-8_dp)
        p(1)  = log(max(omega, 1.0e-12_dp))
        p(2)  = log(max(alpha, 1.0e-12_dp) / slack)
        p(3)  = log(max(beta,  1.0e-12_dp) / slack)
        p(4)  = log(rr / (1.0_dp - rr))
        p(5)  = log(pp / (1.0_dp - pp))
        p(6)  = theta
    end subroutine cnagarch_inv_transform

    ! -------------------------------------------------------------------------
    ! CNAGARCH h/q recursion
    ! -------------------------------------------------------------------------
    subroutine cnagarch_hq(y, params, h_out, q_out)
        real(dp),             intent(in)  :: y(:)
        type(garch_params_t), intent(in)  :: params
        real(dp),             intent(out) :: h_out(:), q_out(:)
        real(dp) :: backcast, h_prev, q_prev, h, qq, y_prev, asym
        integer  :: t

        backcast = max(sum(y**2) / real(size(y), dp), 1.0e-12_dp)
        q_prev   = params%omega / max(1.0_dp - params%extra1, 1.0e-8_dp)
        q_prev   = max(q_prev, 1.0e-12_dp)
        h_prev   = max(q_prev, backcast)
        y_prev   = 0.0_dp
        do t = 1, size(y)
            asym     = (y_prev - params%theta * sqrt(h_prev))**2
            qq       = params%omega + params%extra1*q_prev + params%extra2*(asym - h_prev)
            qq       = max(qq, 1.0e-12_dp)
            h        = qq + params%alpha*(asym - q_prev) + params%beta*(h_prev - q_prev)
            h        = max(h, 1.0e-12_dp)
            h_out(t) = h;  q_out(t) = qq
            q_prev = qq;   h_prev = h;   y_prev = y(t)
        end do
    end subroutine cnagarch_hq

    subroutine cnagarch_skew_kurt(y, params, skew, ekurt)
        real(dp),             intent(in)  :: y(:)
        type(garch_params_t), intent(in)  :: params
        real(dp),             intent(out) :: skew, ekurt
        real(dp) :: hh(size(y)), qq(size(y)), m3, m4, sig2
        integer  :: t, n

        n = size(y)
        call cnagarch_hq(y, params, hh, qq)
        m3 = 0.0_dp;  m4 = 0.0_dp;  sig2 = 0.0_dp
        do t = 1, n
            associate(z => y(t) / sqrt(hh(t)))
                sig2 = sig2 + z**2
                m3   = m3   + z**3
                m4   = m4   + z**4
            end associate
        end do
        sig2 = sig2 / n;  m3 = m3 / n;  m4 = m4 / n
        skew  = m3 / sig2**1.5_dp
        ekurt = m4 / sig2**2 - 3.0_dp
    end subroutine cnagarch_skew_kurt

    ! -------------------------------------------------------------------------
    ! CNAGARCH NLL/n objective (uses cnagarch_data via host association)
    ! -------------------------------------------------------------------------
    real(dp) function cnagarch_value(p)
        real(dp), intent(in) :: p(cnagarch_np)
        type(garch_params_t) :: par
        real(dp) :: hh(size(cnagarch_data)), qq(size(cnagarch_data))
        integer  :: t, n

        call cnagarch_transform(p, par%omega, par%alpha, par%beta, &
                                par%extra1, par%extra2, par%theta)
        n = size(cnagarch_data)
        call cnagarch_hq(cnagarch_data, par, hh, qq)
        cnagarch_value = 0.0_dp
        do t = 1, n
            cnagarch_value = cnagarch_value + log(hh(t)) + cnagarch_data(t)**2 / hh(t)
        end do
        cnagarch_value = 0.5_dp * cnagarch_value / real(n, dp) + log_sqrt_2pi
    end function cnagarch_value

    subroutine cnagarch_obj(p, np_in, f, g)
        integer,  intent(in)  :: np_in
        real(dp), intent(in)  :: p(np_in)
        real(dp), intent(out) :: f, g(np_in)
        real(dp) :: pp(cnagarch_np), fp, fm, step
        integer  :: i

        f = cnagarch_value(p)
        do i = 1, np_in
            step    = 1.0e-5_dp * max(abs(p(i)), 1.0_dp)
            pp      = p;  pp(i) = p(i) + step;  fp = cnagarch_value(pp)
            pp      = p;  pp(i) = p(i) - step;  fm = cnagarch_value(pp)
            g(i)    = (fp - fm) / (2.0_dp * step)
        end do
    end subroutine cnagarch_obj

    ! -------------------------------------------------------------------------
    ! Fit CNAGARCH via BFGS, warm-starting from CSGARCH params (theta=0) and
    ! also trying theta = 0.5, 1.0, -0.5 from the same CSGARCH base.
    ! -------------------------------------------------------------------------
    subroutine fit_cnagarch_local(cs_params, f_best, params_out, niter_best, conv_best)
        type(garch_params_t), intent(in)  :: cs_params
        real(dp),             intent(out) :: f_best
        type(garch_params_t), intent(out) :: params_out
        integer,              intent(out) :: niter_best
        logical,              intent(out) :: conv_best

        real(dp), parameter :: theta_starts(*) = [0.0_dp, 0.5_dp, 1.0_dp, -0.5_dp]
        real(dp) :: p(cnagarch_np), p_best(cnagarch_np), f_try
        integer  :: is, niter_try
        logical  :: conv_try

        f_best   = huge(1.0_dp)
        p_best   = 0.0_dp
        niter_best = 0
        conv_best  = .false.

        do is = 1, size(theta_starts)
            call cnagarch_inv_transform(cs_params%omega, cs_params%alpha, cs_params%beta, &
                                        cs_params%extra1, cs_params%extra2, theta_starts(is), p)
            call bfgs_minimize(cnagarch_obj, p, cnagarch_np, max_iter, gtol, &
                               f_try, niter_try, conv_try)
            if (f_try < f_best) then
                f_best = f_try;  p_best = p;  niter_best = niter_try;  conv_best = conv_try
            end if
        end do

        params_out = garch_params_t()
        call cnagarch_transform(p_best, params_out%omega, params_out%alpha, params_out%beta, &
                                params_out%extra1, params_out%extra2, params_out%theta)
    end subroutine fit_cnagarch_local

end program xcgarch

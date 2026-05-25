! Compare standard NAGARCH against the one-sided news-impact variant on the
! ETF price panel, using demeaned log returns and the best of several starts.
program xcompare_nagarch_news
    use kind_mod,    only: dp
    use csv_mod,     only: read_price_csv, print_price_sample_info
    use nagarch_mod, only: nagarch_set_data, nagarch_set_news_impact, nagarch_obj, &
                           nagarch_transform, nagarch_inv_transform
    use bfgs_mod,    only: bfgs_minimize
    implicit none

    character(len=*), parameter :: prices_file = "spy_efa_eem_tlt_lqd.csv"
    integer, parameter :: max_iter = 1000
    real(dp), parameter :: gtol = 1.0e-7_dp
    integer, parameter :: nstart = 8

    integer, allocatable :: dates(:)
    character(len=32), allocatable :: col_names(:)
    real(dp), allocatable :: prices(:,:), ret(:)
    real(dp) :: starts(4,nstart)
    real(dp) :: p0(4), p(4), f_best
    real(dp) :: omega, alpha, beta, theta
    real(dp) :: logl_std, logl_one, total_delta, mu, theta_std, theta_one
    integer :: nprices, ncols, nobs, icol, best_iter
    logical :: best_conv

    starts(:,1) = [1.0e-5_dp, 0.04_dp, 0.90_dp, -0.50_dp]
    starts(:,2) = [1.0e-5_dp, 0.06_dp, 0.88_dp,  0.00_dp]
    starts(:,3) = [1.0e-5_dp, 0.06_dp, 0.88_dp,  0.50_dp]
    starts(:,4) = [1.0e-5_dp, 0.04_dp, 0.90_dp,  1.00_dp]
    starts(:,5) = [5.0e-6_dp, 0.10_dp, 0.80_dp, -0.50_dp]
    starts(:,6) = [5.0e-6_dp, 0.10_dp, 0.80_dp,  0.00_dp]
    starts(:,7) = [5.0e-6_dp, 0.10_dp, 0.80_dp,  0.50_dp]
    starts(:,8) = [5.0e-6_dp, 0.10_dp, 0.80_dp,  1.00_dp]

    call read_price_csv(prices_file, dates, col_names, prices)
    nprices = size(prices, 1)
    ncols   = size(prices, 2)
    nobs    = nprices - 1
    allocate(ret(nobs))

    call print_price_sample_info(prices_file, dates, ncols)
    write(*,'(A)') "Asset       logL/n std   logL/n one   delta/n    delta total   theta std  theta one"
    write(*,'(A)') repeat("-", 88)

    total_delta = 0.0_dp
    do icol = 1, ncols
        ret = log(prices(2:nprices,icol) / prices(1:nprices-1,icol))
        mu  = sum(ret) / real(nobs, dp)
        ret = ret - mu
        call nagarch_set_data(ret, nobs)

        call fit_mode(.false., f_best, omega, alpha, beta, theta, best_iter, best_conv)
        logl_std = -f_best
        theta_std = theta
        call fit_mode(.true., f_best, omega, alpha, beta, theta, best_iter, best_conv)
        logl_one = -f_best
        theta_one = theta

        total_delta = total_delta + real(nobs, dp) * (logl_one - logl_std)
        write(*,'(A10,2F12.6,F10.6,F14.3,2F11.4)') trim(col_names(icol)), logl_std, logl_one, &
            logl_one - logl_std, real(nobs, dp) * (logl_one - logl_std), &
            theta_std, theta_one
    end do

    write(*,'(A)') repeat("-", 88)
    write(*,'(A,F14.3)') "Total log-likelihood delta, one-sided minus standard: ", total_delta

contains

    ! Fit one NAGARCH news-impact mode from all configured starts and return
    ! the best objective value and structural parameters.
    subroutine fit_mode(zero_above_shift, f_out, omega_out, alpha_out, beta_out, theta_out, iter_out, conv_out)
        logical, intent(in) :: zero_above_shift
        real(dp), intent(out) :: f_out, omega_out, alpha_out, beta_out, theta_out
        integer, intent(out) :: iter_out
        logical, intent(out) :: conv_out
        real(dp) :: f_try
        integer :: istart, iter_try
        logical :: conv_try

        call nagarch_set_news_impact(zero_above_shift)
        f_out = huge(1.0_dp)
        omega_out = 0.0_dp; alpha_out = 0.0_dp; beta_out = 0.0_dp; theta_out = 0.0_dp
        iter_out = 0; conv_out = .false.
        do istart = 1, nstart
            call nagarch_inv_transform(starts(1,istart), starts(2,istart), starts(3,istart), &
                                       starts(4,istart), p0)
            p = p0
            call bfgs_minimize(nagarch_obj, p, 4, max_iter, gtol, f_try, iter_try, conv_try)
            if (f_try < f_out) then
                f_out = f_try
                call nagarch_transform(p, omega_out, alpha_out, beta_out, theta_out)
                iter_out = iter_try
                conv_out = conv_try
            end if
        end do
    end subroutine fit_mode

end program xcompare_nagarch_news

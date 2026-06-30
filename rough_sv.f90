module rough_sv_mod
  use random_mod,  only: random_normal
  use stats_mod,   only: simple_linreg
  use linalg_mod,  only: gauss_elim
  implicit none

  private
  public :: dp, simulate_rough_bergomi, fit_rough_bergomi_returns, estimate_rfsv, &
            rfsv_forecast_weights

  integer, parameter :: dp = kind(1.0d0), nstat = 8

contains

  !> Simulate rough Bergomi-style stochastic volatility price and variance paths.
  !!
  !! The volatility driver is approximated by
  !!
  !!   Y_t = sqrt(2H) int_0^t (t-s)^(H-1/2) dW_s
  !!
  !! and the instantaneous variance is
  !!
  !!   v_t = xi0 * exp(eta * Y_t - 0.5 * eta^2 * t^(2H)).
  !!
  !! The return Brownian increment is correlated with the volatility
  !! Brownian increment using
  !!
  !!   dB_t = rho * dW_t + sqrt(1 - rho^2) * dW_perp_t.
  !!
  !! This simple implementation is O(nstep^2) per path.
  subroutine simulate_rough_bergomi(nstep, npath, tmat, s0, xi0, eta, h, rho, &
                                    s_path, v_path)

    integer, intent(in) :: nstep                    ! Number of time steps.
    integer, intent(in) :: npath                    ! Number of Monte Carlo paths.
    real(dp), intent(in) :: tmat                    ! Maturity or time horizon in years.
    real(dp), intent(in) :: s0                      ! Initial asset price.
    real(dp), intent(in) :: xi0                     ! Initial forward variance level.
    real(dp), intent(in) :: eta                     ! Volatility-of-volatility parameter.
    real(dp), intent(in) :: h                       ! Roughness parameter, in (0, 0.5).
    real(dp), intent(in) :: rho                     ! Return-volatility shock correlation.
    real(dp), intent(out) :: s_path(0:nstep,npath)  ! Simulated asset price paths.
    real(dp), intent(out) :: v_path(0:nstep,npath)  ! Simulated variance paths.

    real(dp), allocatable :: dw(:), dw_perp(:), y(:)
    real(dp) :: dt, sqrt_dt, t, kernel, db, drift_adj, vol, s_old, rho_perp
    integer :: i, j, p

    if (nstep <= 0) error stop "nstep must be positive"
    if (npath <= 0) error stop "npath must be positive"
    if (tmat <= 0.0_dp) error stop "tmat must be positive"
    if (xi0 <= 0.0_dp) error stop "xi0 must be positive"
    if (eta < 0.0_dp) error stop "eta must be non-negative"
    if (h <= 0.0_dp .or. h >= 0.5_dp) error stop "h must be in (0, 0.5)"
    if (abs(rho) > 1.0_dp) error stop "abs(rho) must be <= 1"

    dt = tmat / real(nstep, dp)
    sqrt_dt = sqrt(dt)
    rho_perp = sqrt(max(1.0_dp - rho*rho, 0.0_dp))

    allocate(dw(nstep), dw_perp(nstep), y(nstep))

    do p = 1, npath

      do i = 1, nstep
        dw(i) = sqrt_dt * random_normal()
        dw_perp(i) = sqrt_dt * random_normal()
      end do

      y = 0.0_dp

      do i = 1, nstep
        do j = 1, i
          kernel = (real(i - j + 1, dp) * dt)**(h - 0.5_dp)
          y(i) = y(i) + kernel * dw(j)
        end do

        y(i) = sqrt(2.0_dp * h) * y(i)
      end do

      s_path(0,p) = s0
      v_path(0,p) = xi0

      do i = 1, nstep
        t = real(i, dp) * dt
        drift_adj = -0.5_dp * eta*eta * t**(2.0_dp*h)
        v_path(i,p) = xi0 * exp(eta * y(i) + drift_adj)

        db = rho * dw(i) + rho_perp * dw_perp(i)

        s_old = s_path(i-1,p)
        vol = sqrt(max(v_path(i-1,p), 0.0_dp))

        s_path(i,p) = s_old * exp(-0.5_dp * v_path(i-1,p) * dt + vol * db)
      end do

    end do

    deallocate(dw, dw_perp, y)

  end subroutine simulate_rough_bergomi


  !> Fit rough Bergomi-style parameters to a vector of log returns.
  !!
  !! The method is a simple simulated method of moments estimator.
  !! It estimates xi0 from the annualized return variance and then searches
  !! over supplied grids for eta, h, and rho.
  !!
  !! The fitted parameters minimize the distance between observed and
  !! simulated summary statistics:
  !!
  !!   1. annualized return variance
  !!   2. normalized variance of squared returns
  !!   3. acf(abs(r), lag 1)
  !!   4. acf(abs(r), lag 5)
  !!   5. acf(r^2, lag 1)
  !!   6. acf(r^2, lag 5)
  !!   7. corr(r_t, r_{t+1}^2)
  !!   8. roughness proxy from log rolling realized variance
  subroutine fit_rough_bergomi_returns(ret, dt, npath_fit, eta_grid, h_grid, &
                                       rho_grid, xi0_hat, eta_hat, h_hat, &
                                       rho_hat, obj_min)

    real(dp), intent(in) :: ret(:)         ! Observed log returns.
    real(dp), intent(in) :: dt             ! Time step in years, for example 1/252.
    integer, intent(in) :: npath_fit       ! Number of simulated paths per grid point.
    real(dp), intent(in) :: eta_grid(:)    ! Candidate eta values.
    real(dp), intent(in) :: h_grid(:)      ! Candidate h values.
    real(dp), intent(in) :: rho_grid(:)    ! Candidate rho values.
    real(dp), intent(out) :: xi0_hat       ! Estimated initial variance level.
    real(dp), intent(out) :: eta_hat       ! Estimated volatility-of-volatility.
    real(dp), intent(out) :: h_hat         ! Estimated roughness parameter.
    real(dp), intent(out) :: rho_hat       ! Estimated return-volatility correlation.
    real(dp), intent(out) :: obj_min       ! Minimum simulated moment objective.

    real(dp) :: data_stats(nstat), sim_stats(nstat), obj
    integer :: ie, ih, ir, nret

    nret = size(ret)

    if (nret < 30) error stop "ret must contain at least 30 observations"
    if (dt <= 0.0_dp) error stop "dt must be positive"
    if (npath_fit <= 0) error stop "npath_fit must be positive"
    if (size(eta_grid) == 0) error stop "eta_grid must be nonempty"
    if (size(h_grid) == 0) error stop "h_grid must be nonempty"
    if (size(rho_grid) == 0) error stop "rho_grid must be nonempty"

    call compute_return_stats(ret, dt, data_stats)

    xi0_hat = max(data_stats(1), tiny(1.0_dp))
    eta_hat = eta_grid(1)
    h_hat = h_grid(1)
    rho_hat = rho_grid(1)
    obj_min = huge(1.0_dp)

    do ih = 1, size(h_grid)
      if (h_grid(ih) <= 0.0_dp .or. h_grid(ih) >= 0.5_dp) cycle

      do ie = 1, size(eta_grid)
        if (eta_grid(ie) < 0.0_dp) cycle

        do ir = 1, size(rho_grid)
          if (abs(rho_grid(ir)) > 1.0_dp) cycle

          call average_simulated_stats(nret, npath_fit, dt, xi0_hat, &
                                       eta_grid(ie), h_grid(ih), rho_grid(ir), &
                                       sim_stats)

          obj = moment_objective(data_stats, sim_stats)

          if (obj < obj_min) then
            obj_min = obj
            eta_hat = eta_grid(ie)
            h_hat = h_grid(ih)
            rho_hat = rho_grid(ir)
          end if

        end do
      end do
    end do

  end subroutine fit_rough_bergomi_returns


  !> Compute average simulated return statistics for one parameter vector.
  subroutine average_simulated_stats(nstep, npath, dt, xi0, eta, h, rho, stats)

    integer, intent(in) :: nstep          ! Number of simulated return observations.
    integer, intent(in) :: npath          ! Number of simulated paths.
    real(dp), intent(in) :: dt            ! Time step in years.
    real(dp), intent(in) :: xi0           ! Initial variance level.
    real(dp), intent(in) :: eta           ! Volatility-of-volatility parameter.
    real(dp), intent(in) :: h             ! Roughness parameter.
    real(dp), intent(in) :: rho           ! Return-volatility shock correlation.
    real(dp), intent(out) :: stats(nstat) ! Average simulated summary statistics.

    real(dp), allocatable :: s_path(:,:), v_path(:,:), r(:)
    real(dp) :: stat_one(nstat), tmat, s0
    integer :: i, p

    tmat = real(nstep, dp) * dt
    s0 = 100.0_dp
    stats = 0.0_dp

    allocate(s_path(0:nstep,npath), v_path(0:nstep,npath), r(nstep))

    call simulate_rough_bergomi(nstep, npath, tmat, s0, xi0, eta, h, rho, &
                                s_path, v_path)

    do p = 1, npath
      do i = 1, nstep
        r(i) = log(s_path(i,p) / s_path(i-1,p))
      end do

      call compute_return_stats(r, dt, stat_one)
      stats = stats + stat_one
    end do

    stats = stats / real(npath, dp)

    deallocate(s_path, v_path, r)

  end subroutine average_simulated_stats


  !> Compute summary statistics used by the simulated method of moments fit.
  subroutine compute_return_stats(ret, dt, stats)

    real(dp), intent(in) :: ret(:)        ! Log return vector.
    real(dp), intent(in) :: dt            ! Time step in years.
    real(dp), intent(out) :: stats(nstat) ! Summary statistics.

    real(dp), allocatable :: abs_ret(:), sq_ret(:)
    real(dp) :: v
    integer :: n

    n = size(ret)

    allocate(abs_ret(n), sq_ret(n))

    abs_ret = abs(ret)
    sq_ret = ret*ret
    v = sample_var(ret)

    stats(1) = v / dt
    stats(2) = sample_var(sq_ret) / max(v*v, tiny(1.0_dp))
    stats(3) = acf_lag(abs_ret, 1)
    stats(4) = acf_lag(abs_ret, min(5, n - 1))
    stats(5) = acf_lag(sq_ret, 1)
    stats(6) = acf_lag(sq_ret, min(5, n - 1))
    stats(7) = corr_lag_xy(ret, sq_ret, 1)
    stats(8) = rough_h_proxy(ret, dt)

    deallocate(abs_ret, sq_ret)

  end subroutine compute_return_stats


  !> Compute the weighted distance between observed and simulated statistics.
  pure function moment_objective(data_stats, sim_stats) result(obj)

    real(dp), intent(in) :: data_stats(nstat) ! Observed summary statistics.
    real(dp), intent(in) :: sim_stats(nstat)  ! Simulated summary statistics.
    real(dp) :: obj                           ! Weighted squared distance.

    real(dp), parameter :: w(nstat) = [10.0_dp, 1.0_dp, 5.0_dp, 5.0_dp, &
                                       5.0_dp, 5.0_dp, 5.0_dp, 10.0_dp]
    real(dp) :: scale, diff
    integer :: i

    obj = 0.0_dp

    do i = 1, nstat
      if (i <= 2) then
        scale = max(abs(data_stats(i)), 1.0e-8_dp)
      else
        scale = 1.0_dp
      end if

      diff = (sim_stats(i) - data_stats(i)) / scale
      obj = obj + w(i) * diff * diff
    end do

  end function moment_objective


  !> Estimate a roughness proxy from log rolling realized variance.
  pure function rough_h_proxy(ret, dt) result(h_est)

    real(dp), intent(in) :: ret(:) ! Log return vector.
    real(dp), intent(in) :: dt     ! Time step in years.
    real(dp) :: h_est              ! Roughness proxy.

    real(dp), allocatable :: rv(:), x(:), lx(:), ly(:)
    real(dp) :: msd, slope
    integer :: n, win, m, lag_list(4), i, j, k, nlag

    n = size(ret)

    if (n < 30) then
      h_est = 0.25_dp
      return
    end if

    win = max(5, min(20, n / 10))
    m = n - win + 1

    if (m < 20) then
      h_est = 0.25_dp
      return
    end if

    allocate(rv(m), x(m), lx(4), ly(4))

    do i = 1, m
      rv(i) = sum(ret(i:i+win-1)**2) / max(real(win, dp) * dt, tiny(1.0_dp))
      x(i) = log(max(rv(i), tiny(1.0_dp)))
    end do

    lag_list = [1, 2, 5, 10]
    nlag = 0

    do j = 1, 4
      k = lag_list(j)

      if (k >= m) cycle

      msd = 0.0_dp

      do i = 1, m - k
        msd = msd + (x(i+k) - x(i))**2
      end do

      msd = msd / real(m - k, dp)

      if (msd > 0.0_dp) then
        nlag = nlag + 1
        lx(nlag) = log(real(k, dp) * dt)
        ly(nlag) = log(msd)
      end if
    end do

    if (nlag < 2) then
      h_est = 0.25_dp
    else
      slope = ols_slope(lx(1:nlag), ly(1:nlag))
      h_est = max(0.01_dp, min(0.49_dp, 0.5_dp * slope))
    end if

    deallocate(rv, x, lx, ly)

  end function rough_h_proxy


  !> Compute the sample variance of a vector.
  pure function sample_var(x) result(v)

    real(dp), intent(in) :: x(:) ! Input vector.
    real(dp) :: v               ! Sample variance.

    real(dp) :: xbar
    integer :: n

    n = size(x)

    if (n <= 1) then
      v = 0.0_dp
      return
    end if

    xbar = sum(x) / real(n, dp)
    v = sum((x - xbar)**2) / real(n - 1, dp)

  end function sample_var


  !> Compute the lag-k sample autocorrelation of a vector.
  pure function acf_lag(x, lag) result(acf)

    real(dp), intent(in) :: x(:) ! Input vector.
    integer, intent(in) :: lag   ! Autocorrelation lag.
    real(dp) :: acf              ! Sample autocorrelation.

    real(dp) :: xbar, denom, numer
    integer :: n, i

    n = size(x)

    if (lag <= 0 .or. lag >= n) then
      acf = 0.0_dp
      return
    end if

    xbar = sum(x) / real(n, dp)
    denom = sum((x - xbar)**2)

    if (denom <= tiny(1.0_dp)) then
      acf = 0.0_dp
      return
    end if

    numer = 0.0_dp

    do i = 1, n - lag
      numer = numer + (x(i) - xbar) * (x(i+lag) - xbar)
    end do

    acf = numer / denom

  end function acf_lag


  !> Compute the lagged correlation corr(x_t, y_{t+lag}).
  pure function corr_lag_xy(x, y, lag) result(corr)

    real(dp), intent(in) :: x(:) ! First input vector.
    real(dp), intent(in) :: y(:) ! Second input vector.
    integer, intent(in) :: lag   ! Lag applied to the second vector.
    real(dp) :: corr             ! Lagged correlation.

    real(dp) :: xbar, ybar, sx, sy, sxy
    integer :: n, m, i

    n = min(size(x), size(y))

    if (lag <= 0 .or. lag >= n) then
      corr = 0.0_dp
      return
    end if

    m = n - lag

    xbar = sum(x(1:m)) / real(m, dp)
    ybar = sum(y(1+lag:n)) / real(m, dp)

    sx = 0.0_dp
    sy = 0.0_dp
    sxy = 0.0_dp

    do i = 1, m
      sx = sx + (x(i) - xbar)**2
      sy = sy + (y(i+lag) - ybar)**2
      sxy = sxy + (x(i) - xbar) * (y(i+lag) - ybar)
    end do

    if (sx <= tiny(1.0_dp) .or. sy <= tiny(1.0_dp)) then
      corr = 0.0_dp
    else
      corr = sxy / sqrt(sx * sy)
    end if

  end function corr_lag_xy


  !> Compute the ordinary least squares slope in y = a + b*x.
  pure function ols_slope(x, y) result(b)

    real(dp), intent(in) :: x(:) ! Explanatory variable.
    real(dp), intent(in) :: y(:) ! Response variable.
    real(dp) :: b                ! OLS slope.

    real(dp) :: xbar, ybar, sxx, sxy
    integer :: n

    n = min(size(x), size(y))

    if (n <= 1) then
      b = 0.0_dp
      return
    end if

    xbar = sum(x(1:n)) / real(n, dp)
    ybar = sum(y(1:n)) / real(n, dp)

    sxx = sum((x(1:n) - xbar)**2)
    sxy = sum((x(1:n) - xbar) * (y(1:n) - ybar))

    if (sxx <= tiny(1.0_dp)) then
      b = 0.0_dp
    else
      b = sxy / sxx
    end if

  end function ols_slope


  ! Estimate RFSV parameters from log-MSD regression on daily log realized variance.
  ! Model: log(RV_t) = nu + sigma * W^H(t); slope of log-MSD vs log(lag) gives 2H.
  subroutine estimate_rfsv(x, lags, h_hat, sigma_hat, nu_hat)
    real(dp), intent(in)  :: x(:)        ! daily log realized variance series
    integer,  intent(in)  :: lags(:)     ! lags (in days) used in log-MSD regression
    real(dp), intent(out) :: h_hat       ! estimated Hurst roughness index
    real(dp), intent(out) :: sigma_hat   ! estimated vol of log-RV (fBM diffusion coeff)
    real(dp), intent(out) :: nu_hat      ! estimated mean of log-RV

    integer  :: i, k, n, nlag
    real(dp) :: msd, slope, intercept
    real(dp) :: lx(size(lags)), ly(size(lags))

    n      = size(x)
    nu_hat = sum(x) / real(n, dp)
    nlag   = 0

    do i = 1, size(lags)
      k = lags(i)
      if (k >= n) cycle
      msd = sum((x(k+1:n) - x(1:n-k))**2) / real(n - k, dp)
      if (msd <= 0.0_dp) cycle
      nlag     = nlag + 1
      lx(nlag) = log(real(k, dp))
      ly(nlag) = log(msd)
    end do

    if (nlag < 2) then
      h_hat = 0.1_dp; sigma_hat = 1.0_dp; return
    end if

    call simple_linreg(lx(1:nlag), ly(1:nlag), slope, intercept)
    h_hat     = max(0.01_dp, min(0.49_dp, 0.5_dp * slope))
    sigma_hat = exp(0.5_dp * intercept)
  end subroutine estimate_rfsv

  ! Compute optimal h-step forecast weights for the RFSV model.
  ! Forecast: log(RV_{t+h}) = log(RV_t) + sum_i beta(i)*Delta(t-i+1)
  ! where Delta(t) = log(RV_t) - log(RV_{t-1}).
  ! Weights derived from Yule-Walker equations of the fGn increment process;
  ! sigma cancels so only H is needed.
  subroutine rfsv_forecast_weights(h_step, p, hurst, beta)
    integer,  intent(in)  :: h_step   ! forecast horizon (days)
    integer,  intent(in)  :: p        ! number of past increments used
    real(dp), intent(in)  :: hurst    ! Hurst roughness index H
    real(dp), intent(out) :: beta(p)  ! weights on Delta_t, Delta_{t-1}, ..., Delta_{t-p+1}

    integer  :: i, j, k
    real(dp) :: two_h
    real(dp), allocatable :: rho(:), P_mat(:,:), r_vec(:)

    two_h = 2.0_dp * hurst

    ! fGn normalized autocovariance: rho(0)=1, rho(k)=0.5*((k+1)^{2H}-2k^{2H}+(k-1)^{2H})
    allocate(rho(0:p + h_step - 1))
    rho(0) = 1.0_dp
    do k = 1, p + h_step - 1
      rho(k) = 0.5_dp * (real(k+1,dp)**two_h - 2.0_dp*real(k,dp)**two_h + real(k-1,dp)**two_h)
    end do

    ! Toeplitz correlation matrix P and forecast cross-correlation r
    allocate(P_mat(p,p), r_vec(p))
    do i = 1, p
      do j = 1, p
        P_mat(i,j) = rho(abs(i-j))
      end do
      r_vec(i) = 0.0_dp
      do j = 1, h_step
        r_vec(i) = r_vec(i) + rho(j + i - 1)
      end do
    end do

    call gauss_elim(P_mat, r_vec, p, beta)
    deallocate(rho, P_mat, r_vec)
  end subroutine rfsv_forecast_weights

end module rough_sv_mod

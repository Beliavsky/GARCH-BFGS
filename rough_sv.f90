module rough_sv_mod
  use random_mod, only: random_normal
  implicit none

  integer, parameter :: dp = kind(1.0d0)

contains

  !> Simulate rough Bergomi-style stochastic volatility paths.
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
    real(dp), intent(in) :: h                       ! Roughness parameter, usually in (0, 0.5).
    real(dp), intent(in) :: rho                     ! Correlation between return and volatility shocks.
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

      ! Brownian increments driving volatility and independent return noise.
      do i = 1, nstep
        dw(i) = sqrt_dt * random_normal()
        dw_perp(i) = sqrt_dt * random_normal()
      end do

      ! Build the rough Gaussian driver.
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

        ! Lognormal drift adjustment so E[v_t] is approximately xi0.
        drift_adj = -0.5_dp * eta*eta * t**(2.0_dp*h)

        v_path(i,p) = xi0 * exp(eta * y(i) + drift_adj)

        ! Correlated return Brownian increment.
        db = rho * dw(i) + rho_perp * dw_perp(i)

        s_old = s_path(i-1,p)
        vol = sqrt(max(v_path(i-1,p), 0.0_dp))

        s_path(i,p) = s_old * exp(-0.5_dp * v_path(i-1,p) * dt + vol * db)
      end do

    end do

    deallocate(dw, dw_perp, y)

  end subroutine simulate_rough_bergomi

end module rough_sv_mod


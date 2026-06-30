! Simulate returns from the rough Bergomi stochastic volatility model and
! recover parameters using the simulated method of moments (SMM) estimator.
! Repeats nrep independent replications so estimator variance is visible.
!
! Model:
!   v_t  = xi0 * exp(eta*Y_t - 0.5*eta^2*t^{2H})
!   Y_t  = sqrt(2H) * integral_0^t (t-s)^{H-1/2} dW_s   (Riemann-Liouville fBM)
!   dS/S = sqrt(v_{t-}) * dB_t,  dB_t = rho*dW_t + sqrt(1-rho^2)*dWperp_t
!
! Algorithm notes:
!   Simulation: left-endpoint Volterra sum, O(nstep^2) per path, consistent as dt->0.
!     Drift term -0.5*eta^2*t^{2H} = -0.5*eta^2*Var[Y_t] ensures E[v_t] = xi0.
!   Fitting: SMM grid search over (eta, H, rho); xi0 set to sample ann. variance.
!     8 moments: ann. variance, Var(r^2)/Var(r)^2, ACF(|r|) at lags 1 & 5,
!     ACF(r^2) at lags 1 & 5, corr(r_t, r^2_{t+1}), log-MSD roughness proxy.
!   Limitation: no common random numbers across grid points — MC noise in the
!     objective worsens for small npath_fit.  Increase npath_fit for production.
!
! Usage: xrough_sv.exe

program xrough_sv
    use rough_sv_mod
    use date_mod, only: print_program_header
    implicit none

    ! ---- simulation settings ----
    real(dp), parameter :: tmat      = 5.0_dp            ! time horizon (years)
    real(dp), parameter :: dt        = 1.0_dp / 252.0_dp ! daily time step
    integer,  parameter :: nstep     = nint(tmat / dt)   ! 1260 steps (5 yr daily)
    integer,  parameter :: npath_sim = 1                  ! paths for data generation

    ! ---- true parameters ----
    real(dp), parameter :: xi0_true = 0.04_dp   ! initial variance  (ann. vol ~20%)
    real(dp), parameter :: eta_true = 1.50_dp   ! vol-of-vol
    real(dp), parameter :: h_true   = 0.10_dp   ! roughness  (0 < H < 0.5)
    real(dp), parameter :: rho_true = -0.70_dp  ! leverage correlation
    real(dp), parameter :: s0       = 100.0_dp  ! initial price

    ! ---- SMM fitting settings ----
    integer, parameter :: npath_fit = 20   ! simulated paths per grid point (increase for production)
    integer, parameter :: nrep      = 5   ! independent data replications

    ! ---- parameter grids ----
    integer,  parameter :: n_eta = 6, n_h = 5, n_rho = 5
    real(dp), parameter :: eta_grid(n_eta) = &
        [0.50_dp, 1.00_dp, 1.50_dp, 2.00_dp, 2.50_dp, 3.00_dp]
    real(dp), parameter :: h_grid(n_h) = &
        [0.05_dp, 0.10_dp, 0.15_dp, 0.20_dp, 0.25_dp]
    real(dp), parameter :: rho_grid(n_rho) = &
        [-0.90_dp, -0.70_dp, -0.50_dp, -0.30_dp, -0.10_dp]

    real(dp) :: s_path(0:nstep, npath_sim), v_path(0:nstep, npath_sim)
    real(dp) :: ret(nstep)
    real(dp) :: xi0_fit(nrep), eta_fit(nrep), h_fit(nrep), rho_fit(nrep)
    real(dp) :: obj_min, ann_vol_pct
    real(dp) :: mean_xi0, mean_eta, mean_h, mean_rho
    real(dp) :: std_xi0, std_eta, std_h, std_rho
    integer  :: i, irep
    integer  :: clock_start, clock_end, clock_rate
    real(dp) :: elapsed_s

    call print_program_header("xrough_sv.f90")
    call system_clock(clock_start, clock_rate)

    ! ---- experiment description ----
    write(*,'(A,I0,A,F0.1,A)') &
        "nstep=", nstep, "  T=", tmat, " yr  (daily dt, O(nstep^2) simulation)"
    write(*,'(A,I0,A,I0,A,I0)') &
        "Grid: ", n_eta*n_h*n_rho, " points  npath_fit=", npath_fit, "  nrep=", nrep
    write(*,*)
    write(*,'(A)',advance='no') "eta_grid:"
    do i = 1, n_eta
        write(*,'(1X,F4.2)',advance='no') eta_grid(i)
    end do
    write(*,*)
    write(*,'(A)',advance='no') "  h_grid:"
    do i = 1, n_h
        write(*,'(1X,F4.2)',advance='no') h_grid(i)
    end do
    write(*,*)
    write(*,'(A)',advance='no') "rho_grid:"
    do i = 1, n_rho
        write(*,'(1X,F5.2)',advance='no') rho_grid(i)
    end do
    write(*,'(/)')
    write(*,'(A,F8.4,A,F5.2,A)') "xi0_true =", xi0_true, &
        "  (ann. vol =", sqrt(xi0_true)*100.0_dp, "%)"
    write(*,'(A,F8.4)') "eta_true =", eta_true
    write(*,'(A,F8.4)') "  H_true =", h_true
    write(*,'(A,F8.4)') "rho_true =", rho_true
    write(*,*)

    ! ---- per-replication table ----
    write(*,'(A)') repeat("-", 72)
    write(*,'(A3,1X,A7,2(1X,A7),2(1X,A7),2(1X,A6),2(1X,A7),1X,A9)') &
        "Rep", "AnnVol%", "xi0_tr", "xi0_ft", "eta_tr", "eta_ft", &
        "H_tr", "H_ft", "rho_tr", "rho_ft", "obj_min"
    write(*,'(A)') repeat("-", 72)

    do irep = 1, nrep

        call simulate_rough_bergomi(nstep, npath_sim, tmat, s0, &
            xi0_true, eta_true, h_true, rho_true, s_path, v_path)

        do i = 1, nstep
            ret(i) = log(s_path(i,1) / s_path(i-1,1))
        end do
        ann_vol_pct = sqrt(sum(ret**2) / real(nstep, dp) / dt) * 100.0_dp

        call fit_rough_bergomi_returns(ret, dt, npath_fit, &
            eta_grid, h_grid, rho_grid, &
            xi0_fit(irep), eta_fit(irep), h_fit(irep), rho_fit(irep), obj_min)

        write(*,'(I3,1X,F7.2,2(1X,F7.4),2(1X,F7.4),2(1X,F6.4),2(1X,F7.4),1X,F9.5)') &
            irep, ann_vol_pct, &
            xi0_true, xi0_fit(irep), &
            eta_true, eta_fit(irep), &
            h_true, h_fit(irep), &
            rho_true, rho_fit(irep), &
            obj_min

    end do

    write(*,'(A)') repeat("-", 72)

    ! ---- summary across replications ----
    mean_xi0 = sum(xi0_fit) / real(nrep, dp)
    mean_eta = sum(eta_fit) / real(nrep, dp)
    mean_h   = sum(h_fit)   / real(nrep, dp)
    mean_rho = sum(rho_fit) / real(nrep, dp)

    if (nrep > 1) then
        std_xi0 = sqrt(sum((xi0_fit - mean_xi0)**2) / real(nrep-1, dp))
        std_eta = sqrt(sum((eta_fit - mean_eta)**2) / real(nrep-1, dp))
        std_h   = sqrt(sum((h_fit   - mean_h  )**2) / real(nrep-1, dp))
        std_rho = sqrt(sum((rho_fit - mean_rho)**2) / real(nrep-1, dp))
    else
        std_xi0 = 0.0_dp; std_eta = 0.0_dp
        std_h   = 0.0_dp; std_rho = 0.0_dp
    end if

    write(*,*)
    write(*,'(A)') "Fitted parameter summary:"
    write(*,'(A12,4(1X,A10))') "Statistic", "xi0", "eta", "H", "rho"
    write(*,'(A)') repeat("-", 56)
    write(*,'(A12,4(1X,F10.4))') "True",     xi0_true, eta_true, h_true, rho_true
    write(*,'(A12,4(1X,F10.4))') "Mean_fit", mean_xi0, mean_eta, mean_h, mean_rho
    write(*,'(A12,4(1X,F10.4))') "Std_fit",  std_xi0,  std_eta,  std_h,  std_rho
    write(*,'(A12,4(1X,F10.4))') "Bias", &
        mean_xi0-xi0_true, mean_eta-eta_true, mean_h-h_true, mean_rho-rho_true
    write(*,'(A)') repeat("-", 56)
    write(*,*)
    write(*,'(A)') "Notes:"
    write(*,'(A)') "  xi0_fit = sample ann. variance; varies across reps by design."
    write(*,'(A)') "  eta/H/rho snap to grid points; refine grids or use BFGS polish."
    write(*,'(A)') "  Increase npath_fit (e.g. 100-200) to reduce MC noise in objective."

    call system_clock(clock_end)
    elapsed_s = real(clock_end - clock_start, dp) / real(clock_rate, dp)
    write(*,'(/,A,F8.2,A)') "Elapsed:", elapsed_s, " seconds"

end program xrough_sv

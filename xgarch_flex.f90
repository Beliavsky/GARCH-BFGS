! Test program for garch_flex_mod.
!
! Section 1: Gradient checks for three process/distribution combinations,
!            using GARCH+t simulated data.
!   - proc_garch  + dist_t
!   - proc_egarch + dist_t
!   - proc_gjr    + dist_ged
!
! Section 2: Consistency check — flex(proc_garch, dist_t) vs garch_t_obj
!            on the same data and starting point.  NLL/n values should agree
!            to at least 8 significant figures.
!
! Section 3: Novel fit — flex(proc_egarch, dist_t) on the same data.
!            Prints estimated parameters and NLL/n.

program xgarch_flex
    use kind_mod,        only: dp
    use garch_flex_mod,  only: flex_set_data, flex_set_types, flex_np, flex_obj, &
                                proc_garch, proc_egarch, proc_gjr, &
                                dist_t, dist_ged
    use garch_t_mod,  only: garch_t_simulate, garch_t_set_data, garch_t_obj, &
                                garch_t_transform, garch_t_inv_transform
    use gjr_mod,      only: gjr_inv_transform
    use egarch_mod,   only: egarch_transform, egarch_inv_transform
    use bfgs_mod,     only: bfgs_minimize
    implicit none

    integer,  parameter :: nobs     = 2000
    integer,  parameter :: seed_val = 42
    integer,  parameter :: max_iter = 500
    real(dp), parameter :: gtol     = 1.0e-7_dp
    real(dp), parameter :: h_fd     = 1.0e-5_dp

    ! True simulation parameters (GARCH+t)
    real(dp), parameter :: true_omega = 0.05_dp
    real(dp), parameter :: true_alpha = 0.08_dp
    real(dp), parameter :: true_beta  = 0.85_dp
    real(dp), parameter :: true_nu    = 8.0_dp

    real(dp) :: y(nobs)

    ! Working arrays for gradient checks and fitting
    real(dp) :: p_start_gt(4)   ! starting point for GARCH+t
    real(dp) :: p_start_e(5)    ! starting point for EGARCH+t  (4 proc + 1 dist)
    real(dp) :: p_gjr(5)        ! starting point for GJR+GED   (4 proc + 1 dist)

    real(dp) :: p(5), p_pert(5)
    real(dp) :: f, f_plus, f_minus
    real(dp) :: g(5), g_fd(5)
    real(dp) :: rel_err

    real(dp) :: f_flex, f_gt
    real(dp) :: g_flex(4), g_gt(4)

    real(dp) :: omega_e, alpha_e, gamma_e, beta_e, nu_e
    integer  :: i, np_cur, iter
    logical  :: converged

    ! ================================================================
    ! Simulate GARCH+t data
    ! ================================================================
    call garch_t_simulate(true_omega, true_alpha, true_beta, true_nu, &
                          nobs, seed_val, y)

    ! Load data into the flex module and the garch_t module
    call flex_set_data(y, nobs)
    call garch_t_set_data(y, nobs)

    ! ================================================================
    ! Build starting points
    ! ================================================================
    ! GARCH+t starting point
    call garch_t_inv_transform(0.05_dp, 0.07_dp, 0.83_dp, 10.0_dp, p_start_gt)

    ! EGARCH+t: process params from egarch_inv_transform, dist param reused
    call egarch_inv_transform(-0.05_dp, 0.10_dp, -0.05_dp, 0.95_dp, p_start_e(1:4))
    p_start_e(5) = p_start_gt(4)   ! reuse t-distribution param

    ! GJR+GED: process params from gjr_inv_transform, dist param = log(2)
    call gjr_inv_transform(0.05_dp, 0.07_dp, 0.05_dp, 0.82_dp, p_gjr(1:4))
    p_gjr(5) = 0.7_dp   ! log(2) ≈ 0.693, approx GED shape nu=2

    ! ================================================================
    ! Section 1: Gradient checks
    ! ================================================================
    print '(A)', ""
    print '(A)', " ============================================================"
    print '(A)', " Section 1: Gradient checks (central finite differences, h=1e-5)"
    print '(A)', " ============================================================"

    ! ---- 1a: proc_garch + dist_t ----
    call flex_set_types(proc_garch, dist_t)
    np_cur = flex_np()   ! should be 4
    p(1:np_cur) = p_start_gt(1:np_cur)
    call flex_obj(p(1:np_cur), np_cur, f, g(1:np_cur))
    do i = 1, np_cur
        p_pert(1:np_cur) = p(1:np_cur)
        p_pert(i) = p_pert(i) + h_fd
        call flex_obj(p_pert(1:np_cur), np_cur, f_plus, g_fd(1:np_cur))
        p_pert(1:np_cur) = p(1:np_cur)
        p_pert(i) = p_pert(i) - h_fd
        call flex_obj(p_pert(1:np_cur), np_cur, f_minus, g_fd(1:np_cur))
        g_fd(i) = (f_plus - f_minus) / (2.0_dp * h_fd)
    end do
    print '(A)', ""
    print '(A)', " proc_garch + dist_t:"
    print '(A)', "   i   analytical     fin-diff       rel.error"
    do i = 1, np_cur
        rel_err = abs(g(i) - g_fd(i)) / (abs(g_fd(i)) + 1.0e-15_dp)
        print '(I4,2ES15.6,ES13.3)', i, g(i), g_fd(i), rel_err
    end do

    ! ---- 1b: proc_egarch + dist_t ----
    call flex_set_types(proc_egarch, dist_t)
    np_cur = flex_np()   ! should be 5
    p(1:np_cur) = p_start_e(1:np_cur)
    call flex_obj(p(1:np_cur), np_cur, f, g(1:np_cur))
    do i = 1, np_cur
        p_pert(1:np_cur) = p(1:np_cur)
        p_pert(i) = p_pert(i) + h_fd
        call flex_obj(p_pert(1:np_cur), np_cur, f_plus, g_fd(1:np_cur))
        p_pert(1:np_cur) = p(1:np_cur)
        p_pert(i) = p_pert(i) - h_fd
        call flex_obj(p_pert(1:np_cur), np_cur, f_minus, g_fd(1:np_cur))
        g_fd(i) = (f_plus - f_minus) / (2.0_dp * h_fd)
    end do
    print '(A)', ""
    print '(A)', " proc_egarch + dist_t:"
    print '(A)', "   i   analytical     fin-diff       rel.error"
    do i = 1, np_cur
        rel_err = abs(g(i) - g_fd(i)) / (abs(g_fd(i)) + 1.0e-15_dp)
        print '(I4,2ES15.6,ES13.3)', i, g(i), g_fd(i), rel_err
    end do

    ! ---- 1c: proc_gjr + dist_ged ----
    call flex_set_types(proc_gjr, dist_ged)
    np_cur = flex_np()   ! should be 5
    p(1:np_cur) = p_gjr(1:np_cur)
    call flex_obj(p(1:np_cur), np_cur, f, g(1:np_cur))
    do i = 1, np_cur
        p_pert(1:np_cur) = p(1:np_cur)
        p_pert(i) = p_pert(i) + h_fd
        call flex_obj(p_pert(1:np_cur), np_cur, f_plus, g_fd(1:np_cur))
        p_pert(1:np_cur) = p(1:np_cur)
        p_pert(i) = p_pert(i) - h_fd
        call flex_obj(p_pert(1:np_cur), np_cur, f_minus, g_fd(1:np_cur))
        g_fd(i) = (f_plus - f_minus) / (2.0_dp * h_fd)
    end do
    print '(A)', ""
    print '(A)', " proc_gjr + dist_ged:"
    print '(A)', "   i   analytical     fin-diff       rel.error"
    do i = 1, np_cur
        rel_err = abs(g(i) - g_fd(i)) / (abs(g_fd(i)) + 1.0e-15_dp)
        print '(I4,2ES15.6,ES13.3)', i, g(i), g_fd(i), rel_err
    end do

    ! ================================================================
    ! Section 2: Consistency check — flex(GARCH+t) vs garch_t_obj
    ! ================================================================
    print '(A)', ""
    print '(A)', " ============================================================"
    print '(A)', " Section 2: Consistency check — flex(proc_garch,dist_t) vs garch_t_obj"
    print '(A)', " ============================================================"

    call flex_set_types(proc_garch, dist_t)
    np_cur = flex_np()   ! 4

    ! Evaluate flex at the starting point
    p(1:np_cur) = p_start_gt(1:np_cur)
    call flex_obj(p(1:np_cur), np_cur, f_flex, g_flex(1:np_cur))

    ! Evaluate garch_t_obj at the same starting point
    call garch_t_obj(p_start_gt, 4, f_gt, g_gt)

    print '(A)', ""
    print '(A,ES20.12)', "   flex  NLL/n at start: ", f_flex
    print '(A,ES20.12)', "   gt    NLL/n at start: ", f_gt
    print '(A,ES12.3)',  "   Absolute difference:  ", abs(f_flex - f_gt)

    ! Now fit both with BFGS from the same starting point
    p(1:4) = p_start_gt
    call flex_set_types(proc_garch, dist_t)
    call bfgs_minimize(flex_obj, p(1:4), 4, max_iter, gtol, f_flex, iter, converged)
    print '(A)', ""
    print '(A)', "   After BFGS optimisation:"
    print '(A,ES20.12,A,I0,A,L1)', "   flex  NLL/n: ", f_flex, &
          "   iter: ", iter, "   converged: ", converged

    p(1:4) = p_start_gt
    call bfgs_minimize(garch_t_obj, p(1:4), 4, max_iter, gtol, f_gt, iter, converged)
    print '(A,ES20.12,A,I0,A,L1)', "   gt    NLL/n: ", f_gt, &
          "   iter: ", iter, "   converged: ", converged

    print '(A,ES12.3)', "   Absolute difference:  ", abs(f_flex - f_gt)

    ! ================================================================
    ! Section 3: Novel fit — flex(proc_egarch, dist_t)
    ! ================================================================
    print '(A)', ""
    print '(A)', " ============================================================"
    print '(A)', " Section 3: Novel fit — flex(proc_egarch, dist_t)"
    print '(A)', " ============================================================"

    call flex_set_types(proc_egarch, dist_t)
    np_cur = flex_np()   ! 5
    p(1:np_cur) = p_start_e(1:np_cur)
    call bfgs_minimize(flex_obj, p(1:np_cur), np_cur, max_iter, gtol, f, iter, converged)

    ! Transform process parameters to interpretable form
    call egarch_transform(p(1:4), omega_e, alpha_e, gamma_e, beta_e)
    nu_e = 2.0_dp + 98.0_dp / (1.0_dp + exp(-p(5)))

    print '(A)', ""
    print '(A,I0)', "   nobs = ", nobs
    print '(A)', ""
    print '(A)', "   Estimated parameters (EGARCH+t):"
    print '(A)', "              omega      alpha      gamma       beta         nu"
    print '(A,5F11.4)', "   Estimated ", omega_e, alpha_e, gamma_e, beta_e, nu_e
    print '(A)', ""
    print '(A,I0,A,L1)', "   Iterations: ", iter, "   Converged: ", converged
    print '(A,ES12.4)',  "   NLL/n:      ", f

end program xgarch_flex

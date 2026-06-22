! Select the number of normal mixture components using AIC and BIC.

program xmix_ic
    use date_mod, only: print_program_header
    use normal_mixture_em, only: dp, normal_mixture_type, initialize_mixture, &
                                 fit_mixture, display_parameters, &
                                 generate_mixture_data, log_likelihood, &
                                 smart_initialize, plot_ic
    implicit none

    ! ---- tuneable parameters ----
    integer,  parameter :: n_samples  = 10**3    ! number of observations
    integer,  parameter :: K_MAX      = 6        ! maximum components to try
    integer,  parameter :: K_SIM_MIN  = 1        ! smallest true model to simulate
    integer,  parameter :: K_SIM_MAX  = 4        ! largest true model to simulate
    integer,  parameter :: max_iter   = 200
    real(dp), parameter :: tol        = 1.0e-6_dp

    ! ---- variables ----
    type(normal_mixture_type) :: true_mix, est_mix
    real(dp), allocatable :: data(:)
    real(dp) :: loglik, aic, bic
    real(dp) :: aic_values(K_MAX), bic_values(K_MAX)
    real     :: t_start, t_end
    integer  :: k, k_true, best_k_aic, best_k_bic, n_params, seed_size
    integer  :: chosen_aic(K_SIM_MAX), chosen_bic(K_SIM_MAX)
    integer, allocatable :: seed(:)

    ! ---- random seed ----
    call print_program_header("xmix_ic.f90")
    call random_seed(size=seed_size)
    allocate(seed(seed_size))
    seed = 123456789
    call random_seed(put=seed)

    allocate(data(n_samples))

    ! ---- overall header ----
    print '(A)', "=================================================="
    print '(A)', " Mixture Model Selection via AIC and BIC"
    print '(A)', "=================================================="
    print '(A,I0)', " Observations  : ", n_samples
    print '(A,I0,A,I0)', " True k range  : ", K_SIM_MIN, " .. ", K_SIM_MAX
    print '(A,I0)', " Models tried  : k = 1 .. ", K_MAX
    print '(A)', ""

    call cpu_time(t_start)

    ! ================================================================
    do k_true = K_SIM_MIN, K_SIM_MAX
    ! ================================================================

        ! ---- build true mixture for this k_true ----
        call initialize_mixture(true_mix, k_true)
        select case (k_true)
        case (1)
            true_mix%weights = [1.0_dp]
            true_mix%means   = [0.0_dp]
            true_mix%stdevs  = [2.0_dp]
        case (2)
            true_mix%weights = [0.5_dp, 0.5_dp]
            true_mix%means   = [-3.0_dp, 3.0_dp]
            true_mix%stdevs  = [1.0_dp, 1.0_dp]
        case (3)
            true_mix%weights = [0.3_dp, 0.5_dp, 0.2_dp]
            true_mix%means   = [-5.0_dp, 2.0_dp, 8.0_dp]
            true_mix%stdevs  = [1.5_dp, 1.0_dp, 2.0_dp]
        case (4)
            true_mix%weights = [0.25_dp, 0.25_dp, 0.25_dp, 0.25_dp]
            true_mix%means   = [-6.0_dp, -2.0_dp, 2.0_dp, 6.0_dp]
            true_mix%stdevs  = [1.0_dp,  1.0_dp,  1.0_dp, 1.0_dp]
        end select

        ! ---- generate data ----
        call generate_mixture_data(true_mix, data)

        ! ---- per-simulation header ----
        print '(A)', "=================================================="
        print '(A,I0,A)', " True model: k = ", k_true, " components"
        print '(A)', "=================================================="
        call display_parameters(true_mix, "True Mixture Parameters")

        ! ---- fit k-component models and compute AIC and BIC ----
        !  Parameters per model: k means + k stdevs + (k-1) free weights = 3k-1
        print '(A)', "------------------------------------------------------------"
        print '(A4, A14, A10, A12, A14)', "  k ", "  Log-lik", " #params", "     AIC", "           BIC"
        print '(A)', "------------------------------------------------------------"

        do k = 1, K_MAX
            call initialize_mixture(est_mix, k)
            call smart_initialize(est_mix, data)
            call fit_mixture(est_mix, data, max_iter, tol)

            loglik   = log_likelihood(data, est_mix)
            n_params = 3*k - 1
            aic      = -2.0_dp * loglik + 2.0_dp * n_params
            bic      = -2.0_dp * loglik + log(real(n_samples, dp)) * n_params

            aic_values(k) = aic
            bic_values(k) = bic
            print '(I4, F16.2, I8, F14.2, F14.2)', k, loglik, n_params, aic, bic
        end do

        print '(A)', "------------------------------------------------------------"

        ! ---- find best k by AIC and BIC ----
        best_k_aic = 1
        best_k_bic = 1
        do k = 2, K_MAX
            if (aic_values(k) < aic_values(best_k_aic)) best_k_aic = k
            if (bic_values(k) < bic_values(best_k_bic)) best_k_bic = k
        end do

        print '(A)', ""
        print '(A,I0,A,F14.2)', " Best k (AIC) = ", best_k_aic, "   AIC = ", aic_values(best_k_aic)
        print '(A,I0,A,F14.2)', " Best k (BIC) = ", best_k_bic, "   BIC = ", bic_values(best_k_bic)
        print '(A)', ""

        ! ---- text-based plots ----
        call plot_ic(aic_values, K_MAX, best_k_aic, "AIC")
        call plot_ic(bic_values, K_MAX, best_k_bic, "BIC")

        ! ---- refit and display best AIC model ----
        call initialize_mixture(est_mix, best_k_aic)
        call smart_initialize(est_mix, data)
        call fit_mixture(est_mix, data, max_iter, tol)
        call display_parameters(est_mix, "Estimated Parameters (best AIC model)")

        ! ---- refit and display best BIC model (skip if same k) ----
        if (best_k_bic /= best_k_aic) then
            call initialize_mixture(est_mix, best_k_bic)
            call smart_initialize(est_mix, data)
            call fit_mixture(est_mix, data, max_iter, tol)
            call display_parameters(est_mix, "Estimated Parameters (best BIC model)")
        else
            print '(A)', " (AIC and BIC agree on best k ??? parameters shown above)"
            print '(A)', ""
        end if

        chosen_aic(k_true) = best_k_aic
        chosen_bic(k_true) = best_k_bic

    end do  ! k_true

    call cpu_time(t_end)

    print '(A)', " Model selection summary"
    print '(A)', " ------------------------"
    print '(A)', " True k   AIC chose   BIC chose"
    print '(A)', " ------   ---------   ---------"
    do k_true = K_SIM_MIN, K_SIM_MAX
        print '(I7, I12, I12)', k_true, chosen_aic(k_true), chosen_bic(k_true)
    end do
    print '(A)', ""
    print '(A,F8.3,A)', " CPU time: ", t_end - t_start, " seconds"

end program xmix_ic

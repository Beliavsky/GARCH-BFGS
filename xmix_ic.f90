! Module for normal mixture model fitting using EM algorithm
module normal_mixture_em
    implicit none
    private

    integer, parameter, public :: dp = kind(1.0d0)
    real(kind=dp), parameter :: PI = 3.14159265358979323846_dp

    type, public :: normal_mixture_type
        integer :: n_components
        real(kind=dp), allocatable :: weights(:)
        real(kind=dp), allocatable :: means(:)
        real(kind=dp), allocatable :: stdevs(:)
    end type normal_mixture_type

    public :: initialize_mixture, fit_mixture, display_parameters
    public :: generate_mixture_data, normal_pdf, log_likelihood
    public :: smart_initialize, plot_ic

contains

    subroutine initialize_mixture(mixture, n_components, init_weights, init_means, init_stdevs)
        type(normal_mixture_type), intent(inout) :: mixture
        integer, intent(in) :: n_components
        real(kind=dp), intent(in), optional :: init_weights(:)
        real(kind=dp), intent(in), optional :: init_means(:)
        real(kind=dp), intent(in), optional :: init_stdevs(:)

        mixture%n_components = n_components

        if (allocated(mixture%weights)) deallocate(mixture%weights)
        if (allocated(mixture%means))   deallocate(mixture%means)
        if (allocated(mixture%stdevs))  deallocate(mixture%stdevs)

        allocate(mixture%weights(n_components))
        allocate(mixture%means(n_components))
        allocate(mixture%stdevs(n_components))

        if (present(init_weights)) then
            mixture%weights = init_weights
        else
            mixture%weights = 1.0_dp / n_components
        end if

        if (present(init_means)) then
            mixture%means = init_means
        else
            mixture%means = 0.0_dp
        end if

        if (present(init_stdevs)) then
            mixture%stdevs = init_stdevs
        else
            mixture%stdevs = 1.0_dp
        end if
    end subroutine initialize_mixture

    pure elemental function normal_pdf(x, mu, sigma) result(pdf_val)
        real(kind=dp), intent(in) :: x, mu, sigma
        real(kind=dp) :: pdf_val
        pdf_val = exp(-0.5_dp * ((x - mu) / sigma)**2) / (sigma * sqrt(2.0_dp * PI))
    end function normal_pdf

    pure function mixture_pdf(x, mixture) result(pdf_val)
        real(kind=dp), intent(in) :: x
        type(normal_mixture_type), intent(in) :: mixture
        real(kind=dp) :: pdf_val
        integer :: j

        pdf_val = 0.0_dp
        do j = 1, mixture%n_components
            pdf_val = pdf_val + mixture%weights(j) * &
                      normal_pdf(x, mixture%means(j), mixture%stdevs(j))
        end do
    end function mixture_pdf

    pure function log_likelihood(data, mixture) result(loglik)
        real(kind=dp), intent(in) :: data(:)
        type(normal_mixture_type), intent(in) :: mixture
        real(kind=dp) :: loglik
        real(kind=dp) :: prob
        integer :: i, n

        n = size(data)
        loglik = 0.0_dp
        do i = 1, n
            prob = mixture_pdf(data(i), mixture)
            loglik = loglik + log(prob)
        end do
    end function log_likelihood

    subroutine fit_mixture(mixture, data, max_iter, tol)
        type(normal_mixture_type), intent(inout) :: mixture
        real(kind=dp), intent(in) :: data(:)
        integer, intent(in) :: max_iter
        real(kind=dp), intent(in) :: tol

        real(kind=dp), allocatable :: responsibilities(:,:)
        real(kind=dp) :: old_loglik, new_loglik, diff
        integer :: n, iter

        n = size(data)
        allocate(responsibilities(n, mixture%n_components))

        old_loglik = log_likelihood(data, mixture)

        do iter = 1, max_iter
            call e_step(mixture, data, responsibilities)
            call m_step(mixture, data, responsibilities)
            new_loglik = log_likelihood(data, mixture)
            diff = new_loglik - old_loglik
            if (abs(diff) < tol) exit
            old_loglik = new_loglik
        end do
    end subroutine fit_mixture

    subroutine smart_initialize(mixture, data)
        type(normal_mixture_type), intent(inout) :: mixture
        real(kind=dp), intent(in) :: data(:)
        real(kind=dp) :: data_min, data_max, step, data_mean, data_var
        integer :: i, n

        n = size(data)
        data_min = minval(data)
        data_max = maxval(data)
        step = (data_max - data_min) / (mixture%n_components + 1)

        do i = 1, mixture%n_components
            mixture%means(i) = data_min + i * step
        end do

        mixture%weights = 1.0_dp / mixture%n_components
        data_mean = sum(data) / n
        data_var  = sum((data - data_mean)**2) / n
        mixture%stdevs = sqrt(data_var)
    end subroutine smart_initialize

    subroutine e_step(mixture, data, resp)
        type(normal_mixture_type), intent(in) :: mixture
        real(kind=dp), intent(in) :: data(:)
        real(kind=dp), intent(out) :: resp(:,:)
        real(kind=dp) :: norm_const
        integer :: i, j, n

        n = size(data)
        do i = 1, n
            norm_const = 0.0_dp
            do j = 1, mixture%n_components
                resp(i, j) = mixture%weights(j) * &
                             normal_pdf(data(i), mixture%means(j), mixture%stdevs(j))
                norm_const = norm_const + resp(i, j)
            end do
            if (norm_const > 0.0_dp) then
                resp(i,:) = resp(i,:) / norm_const
            else
                resp(i,:) = 1.0_dp / mixture%n_components
            end if
        end do
    end subroutine e_step

    subroutine m_step(mixture, data, resp)
        type(normal_mixture_type), intent(inout) :: mixture
        real(kind=dp), intent(in) :: data(:)
        real(kind=dp), intent(in) :: resp(:,:)
        real(kind=dp) :: sum_resp, weighted_sum_sq
        integer :: j, n

        n = size(data)
        do j = 1, mixture%n_components
            sum_resp = sum(resp(:, j))
            mixture%weights(j) = sum_resp / n
            if (sum_resp > 0.0_dp) then
                mixture%means(j) = sum(resp(:, j) * data) / sum_resp
                weighted_sum_sq   = sum(resp(:, j) * (data - mixture%means(j))**2)
                mixture%stdevs(j) = sqrt(weighted_sum_sq / sum_resp)
                if (mixture%stdevs(j) < 1.0e-6_dp) mixture%stdevs(j) = 1.0e-6_dp
            end if
        end do
        mixture%weights = mixture%weights / sum(mixture%weights)
    end subroutine m_step

    subroutine display_parameters(mixture, title)
        type(normal_mixture_type), intent(in) :: mixture
        character(len=*), intent(in), optional :: title
        integer :: i

        if (present(title)) then
            print '(A)', trim(title)
            print '(A)', repeat('-', len_trim(title))
        end if

        print '(A5,A10,A15,A15)', 'Comp', 'Weight', 'Mean', 'Std Dev'
        print '(A)', repeat('-', 45)
        do i = 1, mixture%n_components
            print '(I5,F10.4,F15.6,F15.6)', i, mixture%weights(i), &
                  mixture%means(i), mixture%stdevs(i)
        end do
        print *, ''
    end subroutine display_parameters

    subroutine generate_mixture_data(mixture, data)
        type(normal_mixture_type), intent(in) :: mixture
        real(kind=dp), intent(out) :: data(:)
        real(kind=dp) :: u, cumsum, x1, x2, z
        integer :: i, j, component

        do i = 1, size(data)
            call random_number(u)
            cumsum    = 0.0_dp
            component = mixture%n_components
            do j = 1, mixture%n_components
                cumsum = cumsum + mixture%weights(j)
                if (u <= cumsum) then
                    component = j
                    exit
                end if
            end do
            call random_number(x1)
            call random_number(x2)
            z = sqrt(-2.0_dp * log(x1)) * cos(2.0_dp * PI * x2)
            data(i) = mixture%means(component) + mixture%stdevs(component) * z
        end do
    end subroutine generate_mixture_data

    ! Text-based bar plot of an information criterion vs number of components
    subroutine plot_ic(ic_values, k_max, best_k, label)
        real(dp),         intent(in) :: ic_values(:)
        integer,          intent(in) :: k_max, best_k
        character(len=*), intent(in) :: label

        integer,  parameter :: PLOT_WIDTH = 50   ! columns for the bar area

        real(dp) :: ic_min, ic_max, ic_range, frac
        integer  :: k, bar_len
        character(len=PLOT_WIDTH) :: bar
        character(len=1) :: marker

        ic_min   = minval(ic_values(1:k_max))
        ic_max   = maxval(ic_values(1:k_max))
        ic_range = ic_max - ic_min
        if (ic_range < 1.0_dp) ic_range = 1.0_dp   ! guard against flat data

        print '(A)', ""
        print '(A)', " " // trim(label) // " vs Number of Components"
        print '(A)', " (shorter bar = lower " // trim(label) // " = better fit;  * = best k)"
        print '(A)', ""

        do k = 1, k_max
            frac    = (ic_values(k) - ic_min) / ic_range   ! 0 = best, 1 = worst
            bar_len = nint(frac * PLOT_WIDTH)
            if (bar_len < 1) bar_len = 1
            bar = repeat('#', bar_len) // repeat(' ', PLOT_WIDTH - bar_len)

            if (k == best_k) then
                marker = '*'
            else
                marker = ' '
            end if

            print '(A1,I2,A3,A,A,F12.2)', marker, k, ' | ', bar(1:bar_len), ' ', ic_values(k)
        end do

        ! x-axis
        print '(A)', "   +-" // repeat('-', PLOT_WIDTH)
        print '(A)', "     min " // trim(label) // &
                     repeat(' ', PLOT_WIDTH - 9 - len_trim(label)) // "max " // trim(label)
        print '(A)', ""
    end subroutine plot_ic

end module normal_mixture_em

! Main program: for each true k in K_SIM_MIN..K_SIM_MAX, simulate data and
! select best fitting k via AIC and BIC over models with 1..K_MAX components.
program xmix_ic
    use normal_mixture_em, only: dp, normal_mixture_type, initialize_mixture, &
                                 fit_mixture, display_parameters, &
                                 generate_mixture_data, log_likelihood, &
                                 smart_initialize, plot_ic
    implicit none

    ! ---- tuneable parameters ----
    integer,  parameter :: n_samples  = 10**3   ! number of observations
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
            print '(A)', " (AIC and BIC agree on best k — parameters shown above)"
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

! Module for normal mixture model fitting using EM algorithm
module normal_mixture_em
    implicit none
    private
    
    ! Define double precision parameter
    integer, parameter, public :: dp = kind(1.0d0)
    
    ! Define constants
    real(kind=dp), parameter :: PI = 3.14159265358979323846_dp
    
    ! Define derived type for mixture model parameters
    type, public :: normal_mixture_type
        integer :: n_components                ! Number of components
        real(kind=dp), allocatable :: weights(:)     ! Mixing proportions
        real(kind=dp), allocatable :: means(:)       ! Component means
        real(kind=dp), allocatable :: stdevs(:)      ! Component standard deviations
    end type normal_mixture_type
    
    ! Public procedures
    public :: initialize_mixture, fit_mixture, display_parameters
    public :: generate_mixture_data, normal_pdf, log_likelihood
    public :: smart_initialize, plot_ic
    
contains
    ! Initialize mixture parameters
    subroutine initialize_mixture(mixture, n_components, init_weights, init_means, init_stdevs)
        type(normal_mixture_type), intent(inout) :: mixture
        integer, intent(in) :: n_components
        real(kind=dp), intent(in), optional :: init_weights(:)
        real(kind=dp), intent(in), optional :: init_means(:)
        real(kind=dp), intent(in), optional :: init_stdevs(:)
        
        mixture%n_components = n_components
        
        if (allocated(mixture%weights)) deallocate(mixture%weights)
        if (allocated(mixture%means)) deallocate(mixture%means)
        if (allocated(mixture%stdevs)) deallocate(mixture%stdevs)
        
        allocate(mixture%weights(n_components))
        allocate(mixture%means(n_components))
        allocate(mixture%stdevs(n_components))
        
        ! Initialize with default or provided values
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
    
    ! Calculate normal PDF (pure elemental function)
    pure elemental function normal_pdf(x, mu, sigma) result(pdf_val)
        real(kind=dp), intent(in) :: x, mu, sigma
        real(kind=dp) :: pdf_val
        
        pdf_val = exp(-0.5_dp * ((x - mu) / sigma)**2) / (sigma * sqrt(2.0_dp * PI))
    end function normal_pdf
    
    ! Calculate mixture PDF for a single data point (pure function)
    pure function mixture_pdf(x, mixture) result(pdf_val)
        real(kind=dp), intent(in) :: x
        type(normal_mixture_type), intent(in) :: mixture
        real(kind=dp) :: pdf_val
        integer :: j
        
        pdf_val = 0.0_dp
        do j = 1, mixture%n_components
            pdf_val = pdf_val + mixture%weights(j) * normal_pdf(x, mixture%means(j), mixture%stdevs(j))
        end do
    end function mixture_pdf
    
    ! Calculate log-likelihood (pure function)
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
    
    ! EM algorithm to fit mixture model
    subroutine fit_mixture(mixture, data, max_iter, tol, loglik_history)
        type(normal_mixture_type), intent(inout) :: mixture
        real(kind=dp), intent(in) :: data(:)
        integer, intent(in) :: max_iter
        real(kind=dp), intent(in) :: tol
        real(kind=dp), allocatable, intent(out), optional :: loglik_history(:)
        
        real(kind=dp), allocatable :: responsibilities(:,:)
        real(kind=dp) :: old_loglik, new_loglik, diff
        integer :: n, iter
        
        n = size(data)
        
        ! Allocate responsibilities matrix
        allocate(responsibilities(n, mixture%n_components))
        
        ! First log-likelihood calculation
        old_loglik = log_likelihood(data, mixture)
        
        ! Allocate loglik_history if requested
        if (present(loglik_history)) then
            allocate(loglik_history(max_iter+1))
            loglik_history(1) = old_loglik
        end if
        
        ! EM iterations
        do iter = 1, max_iter
            ! E-step: Calculate responsibilities
            call e_step(mixture, data, responsibilities)
            
            ! M-step: Update parameters
            call m_step(mixture, data, responsibilities)
            
            ! Calculate new log-likelihood
            new_loglik = log_likelihood(data, mixture)
            
            ! Store log-likelihood if requested
            if (present(loglik_history)) then
                loglik_history(iter+1) = new_loglik
            end if
            
            ! Check convergence
            diff = new_loglik - old_loglik
            if (abs(diff) < tol) exit
            
            old_loglik = new_loglik
        end do
        
        ! Trim loglik_history if needed
        if (present(loglik_history) .and. iter < max_iter) then
            loglik_history = loglik_history(1:iter+1)
        end if
    end subroutine fit_mixture
    
    ! Smart initialization for EM algorithm
    subroutine smart_initialize(mixture, data)
        type(normal_mixture_type), intent(inout) :: mixture
        real(kind=dp), intent(in) :: data(:)
        real(kind=dp) :: data_min, data_max, step, data_mean, data_var
        integer :: i, n
        
        n = size(data)
        
        ! Simple initialization: spread means uniformly across data range
        data_min = minval(data)
        data_max = maxval(data)
        step = (data_max - data_min) / (mixture%n_components + 1)
        
        do i = 1, mixture%n_components
            mixture%means(i) = data_min + i * step
        end do
        
        ! Initialize weights to be equal
        mixture%weights = 1.0_dp / mixture%n_components
        
        ! Initialize stdevs to data standard deviation
        data_mean = sum(data) / n
        data_var = sum((data - data_mean)**2) / n
        mixture%stdevs = sqrt(data_var)
    end subroutine smart_initialize
    
    ! E-step: Calculate responsibilities
    subroutine e_step(mixture, data, resp)
        type(normal_mixture_type), intent(in) :: mixture
        real(kind=dp), intent(in) :: data(:)
        real(kind=dp), intent(out) :: resp(:,:)
        real(kind=dp) :: norm_const
        integer :: i, j, n
        
        n = size(data)
        
        do i = 1, n
            norm_const = 0.0_dp
            
            ! Calculate unnormalized responsibilities
            do j = 1, mixture%n_components
                resp(i, j) = mixture%weights(j) * normal_pdf(data(i), mixture%means(j), mixture%stdevs(j))
                norm_const = norm_const + resp(i, j)
            end do
            
            ! Normalize responsibilities
            if (norm_const > 0.0_dp) then
                resp(i,:) = resp(i,:) / norm_const
            else
                ! Handle numerical issues
                resp(i,:) = 1.0_dp / mixture%n_components
            end if
        end do
    end subroutine e_step
    
    ! M-step: Update parameters
    subroutine m_step(mixture, data, resp)
        type(normal_mixture_type), intent(inout) :: mixture
        real(kind=dp), intent(in) :: data(:)
        real(kind=dp), intent(in) :: resp(:,:)
        real(kind=dp) :: sum_resp, weighted_sum_sq
        integer :: j, n
        
        n = size(data)
        
        ! Update weights, means, and stdevs for each component
        do j = 1, mixture%n_components
            ! Sum of responsibilities for component j
            sum_resp = sum(resp(:, j))
            
            ! Update weight
            mixture%weights(j) = sum_resp / n
            
            ! Update mean
            if (sum_resp > 0.0_dp) then
                mixture%means(j) = sum(resp(:, j) * data) / sum_resp
            end if
            
            ! Update standard deviation
            if (sum_resp > 0.0_dp) then
                weighted_sum_sq = sum(resp(:, j) * (data - mixture%means(j))**2)
                mixture%stdevs(j) = sqrt(weighted_sum_sq / sum_resp)
                
                ! Avoid zero standard deviation
                if (mixture%stdevs(j) < 1.0e-6_dp) then
                    mixture%stdevs(j) = 1.0e-6_dp
                end if
            end if
        end do
        
        ! Normalize weights to ensure they sum to 1
        mixture%weights = mixture%weights / sum(mixture%weights)
    end subroutine m_step
    
    ! Display parameters in a formatted table
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
            print '(I5,F10.4,F15.6,F15.6)', i, mixture%weights(i), mixture%means(i), mixture%stdevs(i)
        end do
        
        print *, ''  ! Empty line
    end subroutine display_parameters
    
    ! Generate data from a mixture of normals
    subroutine generate_mixture_data(mixture, data)
        type(normal_mixture_type), intent(in) :: mixture
        real(kind=dp), intent(out) :: data(:)
        real(kind=dp) :: u, cumsum, x1, x2, z
        integer :: i, j, component
        
        do i = 1, size(data)
            ! Select component based on weights
            call random_number(u)
            cumsum = 0.0_dp
            component = mixture%n_components  ! Default to last component
            
            do j = 1, mixture%n_components
                cumsum = cumsum + mixture%weights(j)
                if (u <= cumsum) then
                    component = j
                    exit
                end if
            end do
            
            ! Generate normal random variable using Box-Muller transform
            call random_number(x1)
            call random_number(x2)
            z = sqrt(-2.0_dp * log(x1)) * cos(2.0_dp * PI * x2)
            
            ! Scale and shift to match component parameters
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

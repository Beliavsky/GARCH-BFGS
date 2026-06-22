! Test the normal mixture EM algorithm.

program test_mixture_em
    use date_mod, only: print_program_header
    use normal_mixture_em, only: dp, normal_mixture_type, initialize_mixture, &
                                 generate_mixture_data, smart_initialize, &
                                 display_parameters, fit_mixture
    implicit none
    
    ! Parameters
    integer, parameter :: n_samples = 10**6
    integer, parameter :: n_components = 3
    integer, parameter :: max_iter = 100
    real(kind=dp), parameter :: tol = 1.0e-6_dp
    
    ! Variables
    type(normal_mixture_type) :: true_mixture, estimated_mixture
    real(kind=dp), allocatable :: data(:), loglik_history(:)
    integer :: i, seed_size
    integer, allocatable :: seed(:)
    
    ! Initialize random seed
    call print_program_header("xmix.f90")
    call random_seed(size=seed_size)
    allocate(seed(seed_size))
    seed = 123456789  ! Fixed seed for reproducibility
    call random_seed(put=seed)
    
    ! Print simulation parameters
    print '(A)', "Simulation Parameters"
    print '(A)', "---------------------"
    print '(A,I10)', "Number of observations:", n_samples
    print '(A,I6)', "Number of components:", n_components
    print '(A,I6)', "Maximum EM iterations:", max_iter
    print '(A,E10.2)', "Convergence tolerance:", tol
    print *, ""
    
    ! Initialize true mixture model
    call initialize_mixture(true_mixture, n_components)
    true_mixture%weights = [0.3_dp, 0.5_dp, 0.2_dp]
    true_mixture%means = [-5.0_dp, 2.0_dp, 8.0_dp]
    true_mixture%stdevs = [1.5_dp, 1.0_dp, 2.0_dp]
    
    ! Generate data from true mixture
    allocate(data(n_samples))
    call generate_mixture_data(true_mixture, data)
    
    ! Initialize estimated mixture model
    call initialize_mixture(estimated_mixture, n_components)
    
    ! Apply smart initialization before fitting
    call smart_initialize(estimated_mixture, data)
    
    ! Display initial guess before fitting
    call display_parameters(estimated_mixture, "Initial Parameter Guess")
    
    ! Fit the mixture model
    call fit_mixture(estimated_mixture, data, max_iter, tol, loglik_history)
    
    ! Display true and estimated parameters
    call display_parameters(true_mixture, "True Mixture Parameters")
    call display_parameters(estimated_mixture, "Estimated Mixture Parameters")
    
    ! Print log-likelihood history
    print '(A)', "Log-likelihood history (first 10 iterations):"
    do i = 1, min(10, size(loglik_history))
        print '(I5, 1X, F15.6)', i-1, loglik_history(i)  ! Use 1X for explicit space
    end do
    
end program test_mixture_em

module garch_types_mod
    use kind_mod, only: dp
    implicit none
    private

    public :: garch_params_t, garch_fit_result_t

    type :: garch_params_t
        real(dp) :: omega = 0.0_dp
        real(dp) :: alpha = 0.0_dp
        real(dp) :: gamma = 0.0_dp
        real(dp) :: beta  = 0.0_dp
        real(dp) :: theta = 0.0_dp
        real(dp) :: twist = 0.0_dp
        real(dp) :: scale = 0.0_dp
        real(dp) :: extra1 = 0.0_dp
        real(dp) :: extra2 = 0.0_dp
        real(dp), allocatable :: ar_coefs(:)    ! AR(p) mean-equation coefficients phi_1..phi_p
        real(dp), allocatable :: alpha_lags(:)
        real(dp), allocatable :: beta_lags(:)
    end type garch_params_t

    type :: garch_fit_result_t
        character(len=16) :: model = ""
        type(garch_params_t) :: params
        integer :: nparam = 0
        integer :: niter = 0
        logical :: converged = .false.
        real(dp) :: persist = 0.0_dp
        real(dp) :: vol_ann = 0.0_dp
        real(dp) :: logl = 0.0_dp
        real(dp) :: aic = 0.0_dp
        real(dp) :: bic = 0.0_dp
        real(dp) :: skew = 0.0_dp
        real(dp) :: ekurt = 0.0_dp
    end type garch_fit_result_t

end module garch_types_mod
